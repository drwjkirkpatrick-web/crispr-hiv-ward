/*
 * crispr_scan.c — pack FASTA genomes into a compact binary format (.cgb),
 * scan them for CRISPR-Cas9 off-target sites, and map exact k-mer sharing.
 *
 * Build:  gcc -O3 -march=native -pthread -Wall -Wextra -std=c11 crispr_scan.c -o crispr_scan
 *
 * ===========================================================================
 * .CGB PACKED FORMAT (little-endian)
 * ===========================================================================
 *   u32 magic = 0x43474231 ("CGB1")
 *   u32 n_chrom
 *   n_chrom records, in file order:
 *       u32 name_len; name_len bytes (no NUL); u64 chrom_len
 *   then per-chrom data blocks, same order as the records:
 *       ceil(chrom_len/2) bytes, 4 bits per base.
 *   Byte j holds the base at position 2j in the LOW nibble and the base at
 *   position 2j+1 in the HIGH nibble.  Base codes: A=0 C=1 G=2 T=3, anything
 *   else (N, gaps, IUPAC) = 4.  If chrom_len is odd the high nibble of the
 *   last byte is padding (0).  Code 4 "never matches anything": it always
 *   counts as a mismatch inside a window and never satisfies a PAM position
 *   that requires a specific base or base set.
 *
 * ===========================================================================
 * WINDOW/PAM CONVENTIONS (derivation from first principles)
 * ===========================================================================
 * The ordered guide RNA sequence equals the protospacer sequence as read
 * 5'->3' on the PAM-bearing strand, immediately 5' of the PAM.  All window
 * positions below are 0-based coordinates on the plus (reference) strand of
 * the packed chromosome; L is the protospacer length (SPCAS9: 20, SACAS9: 21).
 *
 * PLUS-strand site, window at 0-based p (occupies plus[p .. p+L-1]):
 *   plus[p .. p+L-1] == guide (as given).
 *   The PAM follows on the plus strand at plus[p+L .. p+L+plen-1] and reads
 *   left-to-right (plen = 3 for SPCAS9, 6 for SACAS9).
 *     SPCAS9 NGG:  plus[p+L+1]=='G' && plus[p+L+2]=='G'   (N = plus[p+L], any base)
 *     SPCAS9 NAG:  plus[p+L+1]=='A' && plus[p+L+2]=='G'
 *     SACAS9 NNGRRT: plus[p+L+2]=='G', plus[p+L+3] in {A,G},
 *                    plus[p+L+4] in {A,G}, plus[p+L+5]=='T'
 *     SACAS9 NNGRRV: same first three tests, plus[p+L+5] in {A,C,G}
 *   Reported pam = the actual PAM bases as they read on the plus strand.
 *
 * MINUS-strand site, window at 0-based p (occupies plus[p .. p+L-1]):
 *   plus[p .. p+L-1] == revcomp(guide)  (the protospacer reads revcomp on the
 *   plus strand).  The minus strand runs 5'->3' towards lower plus
 *   coordinates, so the PAM — which is 3' of the protospacer on the
 *   PAM-bearing (minus) strand — occupies plus[p-1 .. p-plen] and reads:
 *       pam[i] = comp(plus[p-1-i]),  i = 0 .. plen-1
 *   Substituting comp() (A<->T, C<->G):
 *     SPCAS9 NGG: comp(plus[p-2])=='G' && comp(plus[p-3])=='G'
 *               i.e. plus[p-2]=='C' && plus[p-3]=='C'; the N of NGG is
 *               comp(plus[p-1]) (any base).
 *     SPCAS9 NAG: plus[p-2]=='T' && plus[p-3]=='C'; pam = comp(plus[p-1])+"AG".
 *     SACAS9 NNGRRT: plus[p-3]=='C', plus[p-4] in {C,T}, plus[p-5] in {C,T},
 *                    plus[p-6]=='A'.
 *     SACAS9 NNGRRV: same first three tests, plus[p-6] in {C,G,T}.
 *   Note NNGRRT and NNGRRV are disjoint at the last base (T vs V={A,C,G});
 *   scanning with both enabled therefore tests the union NNGRR{A,C,G,T},
 *   i.e. any real base at the last position (code 4 still rejected).
 *
 * WORKED EXAMPLE (matches tests/test_scanner.py):
 *   g = AAACCCGGGTTTACGTACGT (SPCAS9, L=20).  A plus-strand NGG site appears
 *   on the plus strand as g + "AGG" (reported pam "AGG").  A minus-strand
 *   NGG site appears on the plus strand as "CC" + "A" + revcomp(g); the PAM
 *   as read on the minus strand is comp('A')+"GG" = "TGG".
 *
 * CODE-4 (N) RULES USED HERE (derived, documented on purpose):
 *   - Window: code 4 vs a guide base is always a mismatch.
 *   - PAM positions requiring a specific base or base set (G, R, T, C, V,
 *     and the SACAS9 union last base) reject code 4.
 *   - PAM pattern positions that accept ANY base (the N of NGG/NAG and the
 *     NN of NNGRRT) accept code 4 as well: the PAM is valid regardless of
 *     which base is actually there, so rejecting it would hide true sites.
 *     Such a base is rendered 'N' in the reported pam string.
 *
 * A site is reported only if the window AND the full PAM lie inside
 * [0, chrom_len).  PAM must match exactly (no mismatches allowed in PAM);
 * the window may have at most max_mm Hamming mismatches vs the guide.
 *
 * ===========================================================================
 * SCAN ALGORITHM (performance notes)
 * ===========================================================================
 * The PAM pattern depends only on the genome, not on the guide, so per tile
 * of TILE_BITS window-start positions we build one bitmap per (system,
 * strand) marking PAM-legal sites, once — amortised over all guides.  Each
 * guide then only walks the set bits of the bitmap and compares its packed
 * window (XOR of 4-bit codes; nonzero-nibble counting; early exit once
 * mismatches exceed max_mm).  Chromosomes are processed by a pthread work
 * queue; each worker holds one chromosome block (pread) plus its hit
 * buffer at a time.  All hits are collected, then globally sorted by
 * (guide order, mismatches, chrom order, position, strand) and written —
 * output is byte-identical for any thread count.
 */

#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdarg.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/* ------------------------------------------------------------------ */
/* constants                                                           */
/* ------------------------------------------------------------------ */

#define MAGIC_CGB1 0x43474231u
#define TILE_BITS 32768u        /* window-start positions per tile        */
#define TILE_WORDS (TILE_BITS / 64u)
#define MAX_THREADS 16
#define SYS_SPCAS9 0u
#define SYS_SACAS9 1u
#define SP_LEN 20               /* SPCAS9 protospacer length             */
#define SP_PAMLEN 3
#define SA_LEN 21               /* SACAS9 protospacer length             */
#define SA_PAMLEN 6

static const char CODE_CHARS[5] = { 'A', 'C', 'G', 'T', 'N' };
static const uint8_t COMP_CODE[5] = { 3, 2, 1, 0, 4 };   /* A<->T C<->G N->N */
static uint8_t CHAR_CODE[256];   /* a/A->0 c/C->1 g/G->2 t/T->3 else 4 */

static void char_code_init(void)
{
    memset(CHAR_CODE, 4, sizeof CHAR_CODE);
    CHAR_CODE['A'] = CHAR_CODE['a'] = 0;
    CHAR_CODE['C'] = CHAR_CODE['c'] = 1;
    CHAR_CODE['G'] = CHAR_CODE['g'] = 2;
    CHAR_CODE['T'] = CHAR_CODE['t'] = 3;
}

static inline bool is_fasta_space(unsigned char c)
{
    return c == ' ' || c == '\t' || c == '\r' || c == '\n' ||
           c == '\v' || c == '\f';
}

/* ------------------------------------------------------------------ */
/* error handling / allocation                                         */
/* ------------------------------------------------------------------ */

static void die(const char *fmt, ...)
    __attribute__((noreturn, format(printf, 1, 2)));
static void die(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fputs("crispr_scan: ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
    exit(1);
}

static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p)
        die("out of memory (%zu bytes)", n);
    return p;
}

static void *xcalloc(size_t n, size_t sz)
{
    void *p = calloc(n ? n : 1, sz ? sz : 1);
    if (!p)
        die("out of memory (%zu x %zu bytes)", n, sz);
    return p;
}

static void *xrealloc(void *p, size_t n)
{
    p = realloc(p, n ? n : 1);
    if (!p)
        die("out of memory (%zu bytes)", n);
    return p;
}

static char *xstrndup(const char *s, size_t n)
{
    char *r = xmalloc(n + 1);
    memcpy(r, s, n);
    r[n] = '\0';
    return r;
}

static int parse_int_range(const char *s, const char *what, int lo, int hi)
{
    char *end = NULL;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno != 0 || end == s || *end != '\0' || v < lo || v > hi)
        die("%s must be an integer in [%d, %d] (got '%s')", what, lo, hi, s);
    return (int)v;
}

static int default_threads(void)
{
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    if (n < 1)
        n = 1;
    if (n > MAX_THREADS)
        n = MAX_THREADS;
    return (int)n;
}

/* ------------------------------------------------------------------ */
/* little-endian I/O helpers                                           */
/* ------------------------------------------------------------------ */

static void wr_bytes(FILE *f, const void *p, size_t n)
{
    if (fwrite(p, 1, n, f) != n)
        die("write error: %s", strerror(errno));
}

static void wr_u32(FILE *f, uint32_t v)
{
    unsigned char b[4] = {
        (unsigned char)v, (unsigned char)(v >> 8),
        (unsigned char)(v >> 16), (unsigned char)(v >> 24)
    };
    wr_bytes(f, b, 4);
}

static void wr_u64(FILE *f, uint64_t v)
{
    unsigned char b[8] = {
        (unsigned char)v, (unsigned char)(v >> 8),
        (unsigned char)(v >> 16), (unsigned char)(v >> 24),
        (unsigned char)(v >> 32), (unsigned char)(v >> 40),
        (unsigned char)(v >> 48), (unsigned char)(v >> 56)
    };
    wr_bytes(f, b, 8);
}

static int rd_u32(FILE *f, uint32_t *v)
{
    unsigned char b[4];
    if (fread(b, 1, 4, f) != 4)
        return -1;
    *v = (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
         ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
    return 0;
}

static int rd_u64(FILE *f, uint64_t *v)
{
    unsigned char b[8];
    if (fread(b, 1, 8, f) != 8)
        return -1;
    *v = 0;
    for (int i = 7; i >= 0; i--)
        *v = (*v << 8) | (uint64_t)b[i];
    return 0;
}

static int pread_full(int fd, void *buf, size_t n, uint64_t off)
{
    uint8_t *p = buf;
    while (n > 0) {
        ssize_t r = pread(fd, p, n, (off_t)off);
        if (r <= 0)
            return -1;
        p += r;
        n -= (size_t)r;
        off += (uint64_t)r;
    }
    return 0;
}

static char *read_whole_file(const char *path, size_t *len_out)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open '%s': %s", path, strerror(errno));
    if (fseeko(f, 0, SEEK_END) != 0)
        die("cannot seek '%s': %s", path, strerror(errno));
    off_t sz = ftello(f);
    if (sz < 0)
        die("cannot stat '%s': %s", path, strerror(errno));
    rewind(f);
    char *buf = xmalloc((size_t)sz + 1);
    if (sz > 0 && fread(buf, 1, (size_t)sz, f) != (size_t)sz)
        die("cannot read '%s'", path);
    buf[sz] = '\0';
    fclose(f);
    *len_out = (size_t)sz;
    return buf;
}

/* ------------------------------------------------------------------ */
/* .cgb container                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    char *name;
    uint64_t len;       /* bases                                          */
    uint64_t off;       /* byte offset of the packed block in the file    */
} Chrom;

typedef struct {
    FILE *f;
    int fd;
    uint32_t n;
    Chrom *c;
} Cgb;

static void cgb_open(const char *path, Cgb *g)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open '%s': %s", path, strerror(errno));
    struct stat st;
    if (fstat(fileno(f), &st) != 0)
        die("cannot stat '%s': %s", path, strerror(errno));
    uint64_t fsize = (uint64_t)st.st_size;

    uint32_t magic = 0;
    if (rd_u32(f, &magic))
        die("'%s': truncated header", path);
    if (magic != MAGIC_CGB1)
        die("'%s': bad magic (not a CGB1 .cgb file)", path);
    uint32_t n = 0;
    if (rd_u32(f, &n))
        die("'%s': truncated header", path);
    if (n > (fsize - 8) / 12 + 1)
        die("'%s': truncated or corrupt header (n_chrom=%u exceeds file size)",
            path, n);

    g->c = xcalloc(n ? n : 1, sizeof(Chrom));
    uint64_t hpos = 8;               /* read position within the header  */
    for (uint32_t i = 0; i < n; i++) {
        uint32_t nlen = 0;
        if (hpos + 4 > fsize || rd_u32(f, &nlen))
            die("'%s': truncated header", path);
        hpos += 4;
        if ((uint64_t)nlen > fsize - hpos)
            die("'%s': truncated header (chrom name)", path);
        char *nm = xmalloc((size_t)nlen + 1);
        if (nlen && fread(nm, 1, nlen, f) != nlen)
            die("'%s': truncated header (chrom name)", path);
        nm[nlen] = '\0';
        hpos += nlen;
        uint64_t len = 0;
        if (hpos + 8 > fsize || rd_u64(f, &len))
            die("'%s': truncated header (chrom length)", path);
        hpos += 8;
        g->c[i].name = nm;
        g->c[i].len = len;
        /* remember the block size for the offset pass below */
        g->c[i].off = (len + 1) / 2;    /* temporarily: block size */
    }
    /* blocks start after all records, in record order */
    uint64_t off = hpos;
    for (uint32_t i = 0; i < n; i++) {
        uint64_t bsz = g->c[i].off;    /* still holding block size */
        if (g->c[i].len > 2 * (fsize - off) + 1)
            die("'%s': truncated (chrom '%s' data extends past end of file)",
                path, g->c[i].name);
        g->c[i].off = off;
        off += bsz;
    }
    g->f = f;
    g->fd = fileno(f);
    g->n = n;
}

static void cgb_close(Cgb *g)
{
    if (g->f)
        fclose(g->f);
    for (uint32_t i = 0; i < g->n; i++)
        free(g->c[i].name);
    free(g->c);
    memset(g, 0, sizeof *g);
}

/* ------------------------------------------------------------------ */
/* packed base access and window primitives                            */
/* ------------------------------------------------------------------ */

static inline uint8_t getb(const uint8_t *blk, uint64_t i)
{
    return (uint8_t)((blk[i >> 1] >> ((i & 1) << 2)) & 0x0Fu);
}

/* packed window: 16 bases in a u64 (nibble i = base i), remainder in a u32 */
typedef struct {
    uint64_t head;   /* bases 0..15                                       */
    uint32_t tail;   /* bases 16..L-1 (nibble j = base 16+j); 0 if L<=16  */
} WinPack;

static WinPack pack_win(const uint8_t *codes, int L)
{
    WinPack w = { 0, 0 };
    for (int i = 0; i < 16 && i < L; i++)
        w.head |= (uint64_t)codes[i] << (4 * i);
    for (int i = 16; i < L; i++)
        w.tail |= (uint32_t)codes[i] << (4 * (i - 16));
    return w;
}

/*
 * Extract genome bases p..p+L-1 in WinPack form.
 * Bounds: the caller guarantees p+L-1 <= chrom_len-1, so all byte reads
 * below stay within the packed block (max byte touched is (p+L-1)/2, and
 * for odd p the helper bytes at (p+15)/2..(p+19)/2 are still inside the
 * window for L >= 20).
 */
static inline WinPack genome_win(const uint8_t *blk, uint64_t p, int L)
{
    WinPack w;
    if (!(p & 1)) {
        uint64_t b = p >> 1;
        memcpy(&w.head, blk + b, 8);
        if (L > 16) {
            uint16_t t;
            memcpy(&t, blk + b + 8, 2);
            w.tail = t;
            if (L == SA_LEN)
                w.tail |= (uint32_t)(blk[b + 10] & 0x0Fu) << 16;
        } else {
            w.tail = 0;
        }
    } else {
        uint64_t b = p >> 1;      /* = (p-1)/2 for odd p */
        uint64_t a;
        memcpy(&a, blk + b, 8);
        uint8_t b8 = blk[b + 8];
        w.head = (a >> 4) | ((uint64_t)b8 << 60);
        if (L > 16) {
            uint8_t b9 = blk[b + 9];
            uint8_t b10 = blk[b + 10];
            w.tail = (uint32_t)(b8 >> 4)
                   | ((uint32_t)(b9 & 0x0Fu) << 4)
                   | ((uint32_t)(b9 >> 4) << 8)
                   | ((uint32_t)(b10 & 0x0Fu) << 12);
            if (L == SA_LEN)
                w.tail |= (uint32_t)(b10 >> 4) << 16;
        } else {
            w.tail = 0;
        }
    }
    return w;
}

/* Count nonzero 4-bit fields of x.  Valid when every nibble of x is <= 7
 * (true here: nibble values are XORs of codes 0..4, which are <= 7). */
static inline unsigned nib_nz64(uint64_t x)
{
    x |= x >> 1;
    x |= x >> 2;
    return (unsigned)__builtin_popcountll(x & 0x1111111111111111ULL);
}

static inline unsigned nib_nz32(uint32_t x)
{
    x |= x >> 1;
    x |= x >> 2;
    return (unsigned)__builtin_popcount(x & 0x11111111u);
}

/* Hamming distance between genome window at p and packed guide g, with
 * early exit once mismatches exceed max_mm.  Code-4 genome bases always
 * XOR nonzero against the ACGT guide and count as mismatches. */
static inline int win_mm(const uint8_t *blk, uint64_t p, const WinPack *g,
                         int L, int max_mm)
{
    WinPack w = genome_win(blk, p, L);
    uint64_t x = w.head ^ g->head;
    int mm = 0;
    if (x) {
        mm = (int)nib_nz64(x);
        if (mm > max_mm)
            return mm;
    }
    uint32_t y = w.tail ^ g->tail;
    if (y) {
        mm += (int)nib_nz32(y);
        if (mm > max_mm)
            return mm;
    }
    return mm;
}

/* ------------------------------------------------------------------ */
/* FASTA parsing                                                       */
/* ------------------------------------------------------------------ */

typedef struct {
    char *name;
    uint64_t len;
} SeqRec;

/* Pass 1: collect record names (first whitespace-delimited token after '>')
 * and sequence lengths.  Sequence characters: whitespace is skipped;
 * A/C/G/T (either case) count as bases; anything else counts as a code-4
 * base.  Blank lines are skipped naturally. */
static SeqRec *fasta_index(FILE *f, const char *path, uint64_t *n_out)
{
    char *line = NULL;
    size_t cap = 0;
    ssize_t ll;
    SeqRec *recs = NULL;
    uint64_t n = 0, capr = 0;
    while ((ll = getline(&line, &cap, f)) != -1) {
        if (line[0] == '>') {
            if (n == capr) {
                capr = capr ? capr * 2 : 16;
                recs = xrealloc(recs, capr * sizeof *recs);
            }
            char *s = line + 1;
            while (*s == ' ' || *s == '\t')
                s++;
            char *e = s;
            while (*e && !is_fasta_space((unsigned char)*e))
                e++;
            recs[n].name = xstrndup(s, (size_t)(e - s));
            recs[n].len = 0;
            n++;
        } else {
            for (ssize_t i = 0; i < ll; i++) {
                if (is_fasta_space((unsigned char)line[i]))
                    continue;
                if (n == 0)
                    die("'%s': sequence data before first FASTA header", path);
                recs[n - 1].len++;
            }
        }
    }
    free(line);
    *n_out = n;
    return recs;
}

/* buffered packed-bytes writer used by pack */
typedef struct {
    FILE *o;
    unsigned char buf[1 << 16];
    size_t wn;
    int half;          /* 1 if a low nibble is pending */
    uint8_t halfval;
} PackOut;

static void po_put(PackOut *P, uint8_t code)
{
    if (!P->half) {
        P->halfval = code;
        P->half = 1;
        return;
    }
    P->half = 0;
    if (P->wn == sizeof P->buf) {
        wr_bytes(P->o, P->buf, P->wn);
        P->wn = 0;
    }
    P->buf[P->wn++] = (unsigned char)(P->halfval | ((unsigned)code << 4));
}

static void po_end_rec(PackOut *P)
{
    if (P->half) {
        P->half = 0;
        if (P->wn == sizeof P->buf) {
            wr_bytes(P->o, P->buf, P->wn);
            P->wn = 0;
        }
        P->buf[P->wn++] = P->halfval;   /* high nibble = 0 padding */
    }
}

static void po_flush(PackOut *P)
{
    if (P->wn) {
        wr_bytes(P->o, P->buf, P->wn);
        P->wn = 0;
    }
}

/* Pass 2: stream the FASTA again, emitting packed data blocks in record
 * order.  Verifies the file did not change between passes. */
static void pack_stream(FILE *f, const SeqRec *recs, uint64_t n, FILE *o)
{
    char *line = NULL;
    size_t cap = 0;
    ssize_t ll;
    int64_t cur = -1;
    uint64_t produced = 0;
    PackOut P;
    memset(&P, 0, sizeof P);
    P.o = o;

    while ((ll = getline(&line, &cap, f)) != -1) {
        if (line[0] == '>') {
            if (cur >= 0) {
                if (produced != recs[cur].len)
                    die("input FASTA changed between passes");
                po_end_rec(&P);
            }
            cur++;
            if ((uint64_t)cur >= n)
                die("input FASTA changed between passes (extra header)");
            produced = 0;
            continue;
        }
        if (cur < 0)
            die("sequence data before first FASTA header");
        for (ssize_t i = 0; i < ll; i++) {
            unsigned char ch = (unsigned char)line[i];
            if (is_fasta_space(ch))
                continue;
            po_put(&P, CHAR_CODE[ch]);
            produced++;
            if (produced > recs[cur].len)
                die("input FASTA changed between passes");
        }
    }
    if (cur >= 0) {
        if (produced != recs[cur].len)
            die("input FASTA changed between passes");
        po_end_rec(&P);
    }
    if (n > 0 && cur != (int64_t)n - 1)
        die("input FASTA changed between passes (missing header)");
    po_flush(&P);
    free(line);
}

/* ------------------------------------------------------------------ */
/* mode: pack                                                          */
/* ------------------------------------------------------------------ */

static void do_pack(const char *in, const char *out)
{
    FILE *fin = fopen(in, "rb");
    if (!fin)
        die("cannot open '%s': %s", in, strerror(errno));
    uint64_t n = 0;
    SeqRec *recs = fasta_index(fin, in, &n);
    if (n > 0xFFFFFFFFu)
        die("'%s': too many sequences (%" PRIu64 ")", in, n);
    fclose(fin);

    FILE *fo = fopen(out, "wb");
    if (!fo)
        die("cannot create '%s': %s", out, strerror(errno));
    wr_u32(fo, MAGIC_CGB1);
    wr_u32(fo, (uint32_t)n);
    for (uint64_t i = 0; i < n; i++) {
        size_t nl = strlen(recs[i].name);
        if (nl > 0xFFFFFFFFu)
            die("chrom name too long");
        wr_u32(fo, (uint32_t)nl);
        wr_bytes(fo, recs[i].name, nl);
        wr_u64(fo, recs[i].len);
    }
    fin = fopen(in, "rb");
    if (!fin)
        die("cannot reopen '%s': %s", in, strerror(errno));
    pack_stream(fin, recs, n, fo);
    if (fflush(fo) != 0 || ferror(fo))
        die("write error on '%s'", out);
    fclose(fin);
    fclose(fo);
    for (uint64_t i = 0; i < n; i++)
        free(recs[i].name);
    free(recs);
}

/* ------------------------------------------------------------------ */
/* mode: unpack                                                        */
/* ------------------------------------------------------------------ */

static void do_unpack(const char *in, const char *out)
{
    Cgb g;
    cgb_open(in, &g);
    FILE *o = fopen(out, "wb");
    if (!o)
        die("cannot create '%s': %s", out, strerror(errno));
    for (uint32_t i = 0; i < g.n; i++) {
        const Chrom *ch = &g.c[i];
        uint64_t bsz = (ch->len + 1) / 2;
        uint8_t *blk = xmalloc(bsz ? bsz : 1);
        if (bsz && pread_full(g.fd, blk, bsz, ch->off))
            die("read error on chrom '%s'", ch->name);
        fprintf(o, ">%s\n", ch->name);
        char linebuf[61];
        size_t col = 0;
        for (uint64_t j = 0; j < ch->len; j++) {
            linebuf[col++] = CODE_CHARS[getb(blk, j)];
            if (col == 60) {
                linebuf[col++] = '\n';
                wr_bytes(o, linebuf, col);
                col = 0;
            }
        }
        if (col) {
            linebuf[col++] = '\n';
            wr_bytes(o, linebuf, col);
        }
        free(blk);
    }
    if (fflush(o) != 0 || ferror(o))
        die("write error on '%s'", out);
    fclose(o);
    cgb_close(&g);
}

/* ------------------------------------------------------------------ */
/* mode: stat                                                          */
/* ------------------------------------------------------------------ */

static void do_stat(const char *path)
{
    Cgb g;
    cgb_open(path, &g);
    uint64_t total = 0;
    for (uint32_t i = 0; i < g.n; i++)
        total += g.c[i].len;
    printf("n_chrom\t%u\n", g.n);
    printf("total_bases\t%" PRIu64 "\n", total);
    for (uint32_t i = 0; i < g.n; i++)
        printf("%s\t%" PRIu64 "\n", g.c[i].name, g.c[i].len);
    cgb_close(&g);
}

/* ------------------------------------------------------------------ */
/* guides.tsv parsing                                                  */
/* ------------------------------------------------------------------ */

typedef struct {
    char *id;
    uint8_t sys;      /* SYS_SPCAS9 / SYS_SACAS9 */
    WinPack fwd;      /* packed guide              */
    WinPack rc;       /* packed revcomp(guide)     */
} Guide;

static void parse_guides(const char *path, Guide **out, uint32_t *nout)
{
    size_t flen = 0;
    char *buf = read_whole_file(path, &flen);
    Guide *gs = NULL;
    uint32_t ng = 0, cap = 0;
    uint64_t lno = 0;
    bool first = true;
    char *p = buf, *end = buf + flen;

    while (p < end) {
        char *nl = memchr(p, '\n', (size_t)(end - p));
        char *le = nl ? nl : end;
        lno++;
        if (le > p && le[-1] == '\r')
            le--;
        size_t llen = (size_t)(le - p);

        bool is_header = first && llen > 0 && p[0] == '#';
        first = false;
        if (!is_header && llen > 0) {
            /* need guide_id \t system \t sequence (extra fields ignored) */
            char *t1 = memchr(p, '\t', llen);
            if (!t1)
                die("guides.tsv line %" PRIu64 ": expected 3 tab-separated "
                    "fields (guide_id, system, sequence)", lno);
            char *t2 = memchr(t1 + 1, '\t', (size_t)(le - (t1 + 1)));
            if (!t2)
                die("guides.tsv line %" PRIu64 ": expected 3 tab-separated "
                    "fields (guide_id, system, sequence)", lno);
            char *t3 = memchr(t2 + 1, '\t', (size_t)(le - (t2 + 1)));
            char *seqend = t3 ? t3 : le;

            size_t idlen = (size_t)(t1 - p);
            size_t syslen = (size_t)(t2 - (t1 + 1));
            size_t seqlen = (size_t)(seqend - (t2 + 1));

            uint8_t sys;
            int L;
            if (syslen == 6 && memcmp(t1 + 1, "SPCAS9", 6) == 0) {
                sys = SYS_SPCAS9;
                L = SP_LEN;
            } else if (syslen == 6 && memcmp(t1 + 1, "SACAS9", 6) == 0) {
                sys = SYS_SACAS9;
                L = SA_LEN;
            } else {
                die("guides.tsv line %" PRIu64 ": unknown system '%.*s' "
                    "(expected SPCAS9 or SACAS9)", lno, (int)syslen, t1 + 1);
            }
            if (seqlen != (size_t)L)
                die("guides.tsv line %" PRIu64 ": guide '%.*s': %s guides "
                    "must be exactly %d nt (got %zu)", lno, (int)idlen, p,
                    sys == SYS_SPCAS9 ? "SPCAS9" : "SACAS9", L, seqlen);

            uint8_t codes[SA_LEN], rccodes[SA_LEN];
            for (size_t i = 0; i < seqlen; i++) {
                uint8_t c = CHAR_CODE[(unsigned char)t2[1 + i]];
                if (c > 3)
                    die("guides.tsv line %" PRIu64 ": guide '%.*s': invalid "
                        "character '%c' in sequence (A/C/G/T only)",
                        lno, (int)idlen, p, t2[1 + i]);
                codes[i] = c;
            }
            for (int i = 0; i < L; i++)
                rccodes[i] = COMP_CODE[codes[L - 1 - i]];

            if (ng == cap) {
                cap = cap ? cap * 2 : 64;
                gs = xrealloc(gs, cap * sizeof *gs);
            }
            gs[ng].id = xstrndup(p, idlen);
            gs[ng].sys = sys;
            gs[ng].fwd = pack_win(codes, L);
            gs[ng].rc = pack_win(rccodes, L);
            ng++;
        }
        p = nl ? nl + 1 : end;
    }
    free(buf);
    *out = gs;
    *nout = ng;
}

/* ------------------------------------------------------------------ */
/* mode: scan                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    uint32_t guide;     /* guide index (input order)                       */
    uint32_t chrom;     /* chrom index (.cgb order)                        */
    uint64_t pos;       /* 0-based plus-strand window start                */
    uint8_t strand;     /* '+' or '-'                                      */
    uint8_t mm;
    char pam[SA_PAMLEN + 1];      /* actual PAM, PAM-bearing strand        */
    char window[SA_LEN + 1];      /* actual plus-strand window              */
} Hit;

typedef struct {
    Hit *v;
    size_t n, cap;
} HitVec;

static void hv_push(HitVec *h, const Hit *x)
{
    if (h->n == h->cap) {
        h->cap = h->cap ? h->cap * 2 : 4096;
        h->v = xrealloc(h->v, h->cap * sizeof *h->v);
    }
    h->v[h->n++] = *x;
}

/* Build a Hit and append.  pam is read off the PAM-bearing strand:
 * plus sites read plus[p+L .. p+L+plen-1] left-to-right; minus sites read
 * comp(plus[p-1-i]) for i = 0..plen-1.  window is always the plus strand. */
static void emit_hit(HitVec *hv, uint32_t gid, uint32_t cid, uint64_t p,
                     int mm, char strand, const uint8_t *blk, int L, int plen)
{
    Hit h;
    memset(&h, 0, sizeof h);
    h.guide = gid;
    h.chrom = cid;
    h.pos = p;
    h.strand = (uint8_t)strand;
    h.mm = (uint8_t)mm;
    if (strand == '+') {
        for (int i = 0; i < plen; i++)
            h.pam[i] = CODE_CHARS[getb(blk, p + (uint64_t)L + (uint64_t)i)];
    } else {
        for (int i = 0; i < plen; i++)
            h.pam[i] = CODE_CHARS[COMP_CODE[getb(blk, p - 1 - (uint64_t)i)]];
    }
    for (int i = 0; i < L; i++)
        h.window[i] = CODE_CHARS[getb(blk, p + (uint64_t)i)];
    hv_push(hv, &h);
}

static int hit_cmp(const void *A, const void *B)
{
    const Hit *a = (const Hit *)A;
    const Hit *b = (const Hit *)B;
    if (a->guide != b->guide)
        return a->guide < b->guide ? -1 : 1;
    if (a->mm != b->mm)
        return (int)a->mm - (int)b->mm;
    if (a->chrom != b->chrom)
        return a->chrom < b->chrom ? -1 : 1;
    if (a->pos != b->pos)
        return a->pos < b->pos ? -1 : 1;
    if (a->strand != b->strand)
        return a->strand == (uint8_t)'+' ? -1 : 1;   /* '+' before '-' */
    return 0;
}

/* ----- per-tile PAM bitmaps (built once per tile, shared by guides) ----- */

/* SPCAS9 plus: window at p (len 20), PAM plus[p+20..p+22] = N G/A G */
static void build_sp_plus(const uint8_t *blk, uint64_t len, uint64_t t0,
                          uint64_t t1, bool ngg_only, uint64_t *bm)
{
    memset(bm, 0, TILE_WORDS * sizeof *bm);
    if (len < SP_LEN + 3)
        return;
    uint64_t p_hi = len - (SP_LEN + 3);
    for (uint64_t p = t0; p < t1 && p <= p_hi; p++) {
        uint8_t b1 = getb(blk, p + SP_LEN + 1);
        uint8_t b2 = getb(blk, p + SP_LEN + 2);
        bool ok = (b1 == 2 && b2 == 2);                    /* NGG */
        if (!ok && !ngg_only)
            ok = (b1 == 0 && b2 == 2);                     /* NAG */
        if (ok)
            bm[(p - t0) >> 6] |= 1ULL << ((p - t0) & 63);
    }
}

/* SPCAS9 minus: window at p, PAM on minus reads comp(plus[p-1]) etc. */
static void build_sp_minus(const uint8_t *blk, uint64_t len, uint64_t t0,
                           uint64_t t1, bool ngg_only, uint64_t *bm)
{
    memset(bm, 0, TILE_WORDS * sizeof *bm);
    if (len < SP_LEN + 3)
        return;
    uint64_t p_hi = len - SP_LEN;
    uint64_t p_lo = t0 < 3 ? 3 : t0;
    for (uint64_t p = p_lo; p < t1 && p <= p_hi; p++) {
        uint8_t b2 = getb(blk, p - 2);
        uint8_t b3 = getb(blk, p - 3);
        bool ok = (b2 == 1 && b3 == 1);                    /* CC = NGG */
        if (!ok && !ngg_only)
            ok = (b2 == 3 && b3 == 1);                     /* TC = NAG */
        if (ok)
            bm[(p - t0) >> 6] |= 1ULL << ((p - t0) & 63);
    }
}

/* SACAS9 plus: window at p (len 21), PAM plus[p+21..p+26] = N N G R R T/V */
static void build_sa_plus(const uint8_t *blk, uint64_t len, uint64_t t0,
                          uint64_t t1, bool canon_only, uint64_t *bm)
{
    memset(bm, 0, TILE_WORDS * sizeof *bm);
    if (len < SA_LEN + 6)
        return;
    uint64_t p_hi = len - (SA_LEN + 6);
    for (uint64_t p = t0; p < t1 && p <= p_hi; p++) {
        if (getb(blk, p + SA_LEN + 2) != 2)                /* G */
            continue;
        uint8_t r1 = getb(blk, p + SA_LEN + 3);
        if (r1 != 0 && r1 != 2)                           /* R */
            continue;
        uint8_t r2 = getb(blk, p + SA_LEN + 4);
        if (r2 != 0 && r2 != 2)                           /* R */
            continue;
        uint8_t lb = getb(blk, p + SA_LEN + 5);
        if (canon_only ? (lb != 3) : (lb > 3))             /* T / V */
            continue;
        bm[(p - t0) >> 6] |= 1ULL << ((p - t0) & 63);
    }
}

/* SACAS9 minus: PAM on minus reads comp(plus[p-1])..comp(plus[p-6]) */
static void build_sa_minus(const uint8_t *blk, uint64_t len, uint64_t t0,
                           uint64_t t1, bool canon_only, uint64_t *bm)
{
    memset(bm, 0, TILE_WORDS * sizeof *bm);
    if (len < SA_LEN + 6)
        return;
    uint64_t p_hi = len - SA_LEN;
    uint64_t p_lo = t0 < 6 ? 6 : t0;
    for (uint64_t p = p_lo; p < t1 && p <= p_hi; p++) {
        if (getb(blk, p - 3) != 1)                         /* C */
            continue;
        uint8_t r1 = getb(blk, p - 4);
        if (r1 != 1 && r1 != 3)                           /* {C,T} */
            continue;
        uint8_t r2 = getb(blk, p - 5);
        if (r2 != 1 && r2 != 3)                           /* {C,T} */
            continue;
        uint8_t lb = getb(blk, p - 6);
        if (canon_only ? (lb != 0) : (lb > 3))             /* A / {C,G,T} */
            continue;
        bm[(p - t0) >> 6] |= 1ULL << ((p - t0) & 63);
    }
}

/* Walk the two bitmaps for one guide and one tile. */
static void scan_system_sites(const uint8_t *blk, uint64_t t0,
                              const uint64_t *bmp, const uint64_t *bmm,
                              const Guide *gd, uint32_t gid, uint32_t cid,
                              int max_mm, int L, int plen, HitVec *hv)
{
    for (uint32_t wi = 0; wi < TILE_WORDS; wi++) {
        uint64_t wp = bmp[wi];
        while (wp) {
            unsigned b = (unsigned)__builtin_ctzll(wp);
            wp &= wp - 1;
            uint64_t p = t0 + ((uint64_t)wi << 6) + b;
            int mm = win_mm(blk, p, &gd->fwd, L, max_mm);
            if (mm <= max_mm)
                emit_hit(hv, gid, cid, p, mm, '+', blk, L, plen);
        }
        uint64_t wm = bmm[wi];
        while (wm) {
            unsigned b = (unsigned)__builtin_ctzll(wm);
            wm &= wm - 1;
            uint64_t p = t0 + ((uint64_t)wi << 6) + b;
            int mm = win_mm(blk, p, &gd->rc, L, max_mm);
            if (mm <= max_mm)
                emit_hit(hv, gid, cid, p, mm, '-', blk, L, plen);
        }
    }
}

typedef struct {
    const Cgb *g;
    const Guide *guides;
    uint32_t ng;
    int max_mm;
    bool ngg_only, canon_only;
    bool has_sp, has_sa;
    atomic_uint next;
    HitVec hits[MAX_THREADS];
    uint64_t *bm[MAX_THREADS];
    unsigned spawned;
} ScanShared;

typedef struct {
    ScanShared *sh;
    unsigned tid;
} WorkerArg;

static void scan_chrom(ScanShared *sh, unsigned tid, uint32_t cid)
{
    const Chrom *ch = &sh->g->c[cid];
    uint64_t len = ch->len;
    if (len == 0)
        return;
    uint64_t bsz = (len + 1) / 2;
    uint8_t *blk = xmalloc(bsz);
    if (pread_full(sh->g->fd, blk, bsz, ch->off))
        die("I/O error reading chrom '%s'", ch->name);
    uint64_t *bm = sh->bm[tid];
    HitVec *hv = &sh->hits[tid];

    uint64_t ntiles = (len + TILE_BITS - 1) / TILE_BITS;
    for (uint64_t t = 0; t < ntiles; t++) {
        uint64_t t0 = t * TILE_BITS;
        uint64_t t1 = t0 + TILE_BITS;
        if (t1 > len)
            t1 = len;
        if (sh->has_sp) {
            build_sp_plus(blk, len, t0, t1, sh->ngg_only, bm);
            build_sp_minus(blk, len, t0, t1, sh->ngg_only, bm + TILE_WORDS);
        }
        if (sh->has_sa) {
            build_sa_plus(blk, len, t0, t1, sh->canon_only,
                          bm + 2 * TILE_WORDS);
            build_sa_minus(blk, len, t0, t1, sh->canon_only,
                           bm + 3 * TILE_WORDS);
        }
        for (uint32_t gi = 0; gi < sh->ng; gi++) {
            const Guide *gd = &sh->guides[gi];
            if (gd->sys == SYS_SPCAS9) {
                scan_system_sites(blk, t0, bm, bm + TILE_WORDS, gd, gi, cid,
                                  sh->max_mm, SP_LEN, SP_PAMLEN, hv);
            } else {
                scan_system_sites(blk, t0, bm + 2 * TILE_WORDS,
                                  bm + 3 * TILE_WORDS, gd, gi, cid,
                                  sh->max_mm, SA_LEN, SA_PAMLEN, hv);
            }
        }
    }
    free(blk);
}

static void *scan_worker(void *arg)
{
    WorkerArg *wa = (WorkerArg *)arg;
    ScanShared *sh = wa->sh;
    unsigned tid = wa->tid;
    sh->bm[tid] = xcalloc(4 * TILE_WORDS, sizeof(uint64_t));
    for (;;) {
        uint32_t c = atomic_fetch_add(&sh->next, 1);
        if (c >= sh->g->n)
            break;
        scan_chrom(sh, tid, c);
    }
    free(sh->bm[tid]);
    sh->bm[tid] = NULL;
    return NULL;
}

static const char *sys_name(uint8_t sys)
{
    return sys == SYS_SPCAS9 ? "SPCAS9" : "SACAS9";
}

static void do_scan(int argc, char **argv)
{
    const char *pos[5] = { NULL };
    int npos = 0;
    bool ngg_only = false, canon_only = false;
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "--ngg-only") == 0) {
            ngg_only = true;
            continue;
        }
        if (strcmp(argv[i], "--canonical-only") == 0) {
            canon_only = true;
            continue;
        }
        if (argv[i][0] == '-' && argv[i][1] == '-')
            die("scan: unknown option '%s'", argv[i]);
        if (npos >= 5)
            die("scan: too many arguments (try --help)");
        pos[npos++] = argv[i];
    }
    if (npos < 4)
        die("scan usage: crispr_scan scan <genome.cgb> <guides.tsv> "
            "<max_mm> <out.tsv> [threads]");
    int max_mm = parse_int_range(pos[2], "max_mm", 0, 255);
    int threads = (npos >= 5)
        ? parse_int_range(pos[4], "threads", 1, 1024)
        : default_threads();
    if (threads > MAX_THREADS)
        threads = MAX_THREADS;

    Guide *guides = NULL;
    uint32_t ng = 0;
    parse_guides(pos[1], &guides, &ng);
    bool has_sp = false, has_sa = false;
    for (uint32_t i = 0; i < ng; i++) {
        if (guides[i].sys == SYS_SPCAS9)
            has_sp = true;
        else
            has_sa = true;
    }

    Cgb g;
    cgb_open(pos[0], &g);

    ScanShared sh;
    memset(&sh, 0, sizeof sh);
    sh.g = &g;
    sh.guides = guides;
    sh.ng = ng;
    sh.max_mm = max_mm;
    sh.ngg_only = ngg_only;
    sh.canon_only = canon_only;
    sh.has_sp = has_sp;
    sh.has_sa = has_sa;
    atomic_init(&sh.next, 0);

    unsigned spawn = 0;
    if (ng > 0 && g.n > 0) {
        spawn = (unsigned)threads < g.n ? (unsigned)threads : g.n;
    }
    sh.spawned = spawn;
    pthread_t th[MAX_THREADS];
    WorkerArg wa[MAX_THREADS];
    for (unsigned i = 0; i < spawn; i++) {
        wa[i].sh = &sh;
        wa[i].tid = i;
        if (pthread_create(&th[i], NULL, scan_worker, &wa[i]) != 0)
            die("cannot create thread");
    }
    for (unsigned i = 0; i < spawn; i++)
        pthread_join(th[i], NULL);

    /* merge, sort, write */
    uint64_t total = 0;
    for (unsigned i = 0; i < sh.spawned; i++)
        total += sh.hits[i].n;
    Hit *all = NULL;
    if (total) {
        all = xmalloc((size_t)total * sizeof *all);
        uint64_t k = 0;
        for (unsigned i = 0; i < sh.spawned; i++) {
            memcpy(all + k, sh.hits[i].v, sh.hits[i].n * sizeof *all);
            k += sh.hits[i].n;
            free(sh.hits[i].v);
        }
        qsort(all, (size_t)total, sizeof *all, hit_cmp);
    }

    FILE *o = fopen(pos[3], "wb");
    if (!o)
        die("cannot create '%s': %s", pos[3], strerror(errno));
    fputs("guide_id\tsystem\tchrom\tpos\tstrand\tpam\tmismatches\t"
          "plus_window\n", o);
    for (uint64_t i = 0; i < total; i++) {
        const Hit *h = &all[i];
        fprintf(o, "%s\t%s\t%s\t%" PRIu64 "\t%c\t%s\t%u\t%s\n",
                guides[h->guide].id, sys_name(guides[h->guide].sys),
                g.c[h->chrom].name, h->pos + 1, (char)h->strand, h->pam,
                (unsigned)h->mm, h->window);
    }
    if (fflush(o) != 0 || ferror(o))
        die("write error on '%s'", pos[3]);
    fclose(o);
    free(all);
    for (uint32_t i = 0; i < ng; i++)
        free(guides[i].id);
    free(guides);
    cgb_close(&g);
}

/* ------------------------------------------------------------------ */
/* mode: kmap                                                          */
/* ------------------------------------------------------------------ */

typedef struct {
    uint64_t key;       /* 2-bit packed k-mer code                          */
    uint32_t head;      /* insts index + 1 of first instance, 0 = none      */
    uint8_t used;
} KSlot;

typedef struct {
    uint32_t seq;       /* query sequence index (query order)               */
    uint32_t qpos;      /* 0-based k-mer start in the query sequence        */
    uint8_t orient;     /* 0 = FWD, 1 = RC                                  */
    uint32_t next;      /* insts index + 1 of next instance, 0 = end        */
} KInst;

typedef struct {
    uint32_t seq;
    uint32_t qpos;
    uint8_t orient;
    uint32_t chrom;
    uint64_t cpos;      /* 1-based plus-strand match position               */
    uint64_t code;      /* k-mer as it reads on the genome plus strand      */
} KHit;

typedef struct {
    KHit *v;
    size_t n, cap;
} KHitVec;

static void khv_push(KHitVec *h, const KHit *x)
{
    if (h->n == h->cap) {
        h->cap = h->cap ? h->cap * 2 : 4096;
        h->v = xrealloc(h->v, h->cap * sizeof *h->v);
    }
    h->v[h->n++] = *x;
}

static int khit_cmp(const void *A, const void *B)
{
    const KHit *a = (const KHit *)A;
    const KHit *b = (const KHit *)B;
    if (a->seq != b->seq)
        return a->seq < b->seq ? -1 : 1;
    if (a->qpos != b->qpos)
        return a->qpos < b->qpos ? -1 : 1;
    if (a->orient != b->orient)
        return a->orient < b->orient ? -1 : 1;   /* FWD before RC */
    if (a->chrom != b->chrom)
        return a->chrom < b->chrom ? -1 : 1;
    if (a->cpos != b->cpos)
        return a->cpos < b->cpos ? -1 : 1;
    return 0;
}

typedef struct {
    uint64_t code;      /* 2-bit packed k-mer code                          */
    uint32_t seq;       /* query sequence index (query order)               */
    uint32_t qpos;      /* 0-based k-mer start in the query sequence        */
    uint8_t orient;     /* 0 = FWD, 1 = RC                                  */
} Pending;

/*
 * One streaming pass over the query FASTA.  For every ACGT k-mer at every
 * position, queue two hash entries: the k-mer as given (FWD) and its
 * reverse complement (RC).  k-mers containing non-ACGT characters are
 * skipped (a code-4 base resets both rolling codes).  Returns the record
 * names in file order.
 */
static void kmap_collect(const char *path, int k, uint64_t code_mask,
                         Pending **pend_out, uint64_t *npend_out,
                         char ***names_out, uint64_t *nnames_out)
{
    FILE *f = fopen(path, "rb");
    if (!f)
        die("cannot open '%s': %s", path, strerror(errno));
    char **names = NULL;
    uint64_t nn = 0, ncap = 0;
    Pending *pend = NULL;
    uint64_t np = 0, pcap = 0;
    char *line = NULL;
    size_t lcap = 0;
    ssize_t ll;
    int64_t seq = -1;
    uint64_t ri = 0, fo = 0, rc = 0;
    uint32_t cnt = 0;
    int shift = 2 * (k - 1);

    while ((ll = getline(&line, &lcap, f)) != -1) {
        if (line[0] == '>') {
            if (nn == ncap) {
                ncap = ncap ? ncap * 2 : 16;
                names = xrealloc(names, ncap * sizeof *names);
            }
            char *s = line + 1;
            while (*s == ' ' || *s == '\t')
                s++;
            char *e = s;
            while (*e && !is_fasta_space((unsigned char)*e))
                e++;
            names[nn++] = xstrndup(s, (size_t)(e - s));
            seq++;
            ri = 0;
            fo = 0;
            rc = 0;
            cnt = 0;
            continue;
        }
        if (seq < 0)
            die("'%s': sequence data before first FASTA header", path);
        for (ssize_t i = 0; i < ll; i++) {
            unsigned char ch = (unsigned char)line[i];
            if (is_fasta_space(ch))
                continue;
            uint8_t code = CHAR_CODE[ch];
            if (code > 3) {
                cnt = 0;
                fo = 0;
                rc = 0;
                ri++;
                continue;
            }
            fo = ((fo << 2) | (uint64_t)code) & code_mask;
            rc = (rc >> 2) | ((uint64_t)COMP_CODE[code] << shift);
            cnt++;
            ri++;
            if (cnt >= (uint32_t)k) {
                if (np + 2 > pcap) {
                    pcap = pcap ? pcap * 2 : 4096;
                    pend = xrealloc(pend, pcap * sizeof *pend);
                }
                pend[np].code = fo;
                pend[np].seq = (uint32_t)seq;
                pend[np].qpos = (uint32_t)(ri - (uint64_t)k);
                pend[np].orient = 0;
                np++;
                pend[np].code = rc;
                pend[np].seq = (uint32_t)seq;
                pend[np].qpos = (uint32_t)(ri - (uint64_t)k);
                pend[np].orient = 1;
                np++;
            }
        }
    }
    free(line);
    fclose(f);
    *pend_out = pend;
    *npend_out = np;
    *names_out = names;
    *nnames_out = nn;
}

typedef struct {
    const Cgb *g;
    const KSlot *slots;
    uint64_t slot_mask;
    int slot_shift;
    const KInst *insts;
    uint64_t code_mask;
    int k;
    atomic_uint next;
    KHitVec hits[MAX_THREADS];
    unsigned spawned;
} KmapShared;

typedef struct {
    KmapShared *sh;
    unsigned tid;
} KWorkerArg;

/*
 * Stream one chromosome with a rolling 2-bit k-mer code; a code-4 base
 * resets the window.  Every valid genome k-mer is looked up in the query
 * hash; on a hit every instance under that code emits a row.
 */
static void kmap_chrom(KmapShared *sh, unsigned tid, uint32_t cid)
{
    const Chrom *ch = &sh->g->c[cid];
    uint64_t len = ch->len;
    if (len == 0)
        return;
    uint64_t bsz = (len + 1) / 2;
    uint8_t *blk = xmalloc(bsz);
    if (pread_full(sh->g->fd, blk, bsz, ch->off))
        die("I/O error reading chrom '%s'", ch->name);
    KHitVec *hv = &sh->hits[tid];
    uint64_t f = 0;
    uint32_t cnt = 0;
    for (uint64_t i = 0; i < len; i++) {
        uint8_t b = getb(blk, i);
        if (b > 3) {
            cnt = 0;
            f = 0;
            continue;
        }
        f = ((f << 2) | (uint64_t)b) & sh->code_mask;
        if (++cnt < (uint32_t)sh->k)
            continue;
        uint64_t j = (f * 0x9E3779B97F4A7C15ULL) >> sh->slot_shift;
        for (;;) {
            const KSlot *s = &sh->slots[j & sh->slot_mask];
            if (!s->used)
                break;
            if (s->key == f) {
                uint64_t cpos0 = i + 1 - (uint64_t)sh->k;
                for (uint32_t e = s->head; e; e = sh->insts[e - 1].next) {
                    const KInst *in = &sh->insts[e - 1];
                    KHit h = { in->seq, in->qpos, in->orient, cid,
                               cpos0 + 1, f };
                    khv_push(hv, &h);
                }
                break;
            }
            j++;
        }
    }
    free(blk);
}

static void *kmap_worker(void *arg)
{
    KWorkerArg *wa = (KWorkerArg *)arg;
    KmapShared *sh = wa->sh;
    unsigned tid = wa->tid;
    for (;;) {
        uint32_t c = atomic_fetch_add(&sh->next, 1);
        if (c >= sh->g->n)
            break;
        kmap_chrom(sh, tid, c);
    }
    return NULL;
}

static void do_kmap(int argc, char **argv)
{
    const char *pos[5] = { NULL };
    int npos = 0;
    for (int i = 2; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == '-')
            die("kmap: unknown option '%s'", argv[i]);
        if (npos >= 5)
            die("kmap: too many arguments (try --help)");
        pos[npos++] = argv[i];
    }
    if (npos < 4)
        die("kmap usage: crispr_scan kmap <genome.cgb> <query.fa> <k> "
            "<out.tsv> [threads]");
    int k = parse_int_range(pos[2], "k", 1, 31);
    int threads = (npos >= 5)
        ? parse_int_range(pos[4], "threads", 1, 1024)
        : default_threads();
    if (threads > MAX_THREADS)
        threads = MAX_THREADS;

    uint64_t code_mask = (1ULL << (2 * k)) - 1;    /* k <= 31 -> <= 62 bits */

    Cgb g;
    cgb_open(pos[0], &g);

    Pending *pend = NULL;
    uint64_t npend = 0;
    char **names = NULL;
    uint64_t nnames = 0;
    kmap_collect(pos[1], k, code_mask, &pend, &npend, &names, &nnames);

    /* build the hash table of query k-mers (both orientations) */
    uint64_t cap = 64;
    while (cap < npend * 2)
        cap <<= 1;
    KSlot *slots = xcalloc(cap, sizeof *slots);
    KInst *insts = npend ? xmalloc(npend * sizeof *insts) : NULL;
    uint64_t mask = cap - 1;
    int shift = 64 - __builtin_ctzll(cap);
    for (uint64_t i = 0; i < npend; i++) {
        uint64_t j = (pend[i].code * 0x9E3779B97F4A7C15ULL) >> shift;
        for (;;) {
            KSlot *s = &slots[j & mask];
            if (!s->used) {
                s->used = 1;
                s->key = pend[i].code;
                s->head = (uint32_t)(i + 1);
                insts[i].next = 0;
                break;
            }
            if (s->key == pend[i].code) {
                insts[i].next = s->head;
                s->head = (uint32_t)(i + 1);
                break;
            }
            j++;
        }
        insts[i].seq = pend[i].seq;
        insts[i].qpos = pend[i].qpos;
        insts[i].orient = pend[i].orient;
    }
    free(pend);

    KmapShared sh;
    memset(&sh, 0, sizeof sh);
    sh.g = &g;
    sh.slots = slots;
    sh.slot_mask = mask;
    sh.slot_shift = shift;
    sh.insts = insts;
    sh.code_mask = code_mask;
    sh.k = k;
    atomic_init(&sh.next, 0);

    unsigned spawn = 0;
    if (npend > 0 && g.n > 0) {
        spawn = (unsigned)threads < g.n ? (unsigned)threads : g.n;
    }
    sh.spawned = spawn;
    pthread_t th[MAX_THREADS];
    KWorkerArg ka[MAX_THREADS];
    for (unsigned i = 0; i < spawn; i++) {
        ka[i].sh = &sh;
        ka[i].tid = i;
        if (pthread_create(&th[i], NULL, kmap_worker, &ka[i]) != 0)
            die("cannot create thread");
    }
    for (unsigned i = 0; i < spawn; i++)
        pthread_join(th[i], NULL);

    /* merge, sort, write */
    uint64_t total = 0;
    for (unsigned i = 0; i < sh.spawned; i++)
        total += sh.hits[i].n;
    KHit *all = NULL;
    if (total) {
        all = xmalloc((size_t)total * sizeof *all);
        uint64_t kk = 0;
        for (unsigned i = 0; i < sh.spawned; i++) {
            memcpy(all + kk, sh.hits[i].v, sh.hits[i].n * sizeof *all);
            kk += sh.hits[i].n;
            free(sh.hits[i].v);
        }
        qsort(all, (size_t)total, sizeof *all, khit_cmp);
    }

    FILE *o = fopen(pos[3], "wb");
    if (!o)
        die("cannot create '%s': %s", pos[3], strerror(errno));
    fputs("query_name\tquery_pos\torient\tchrom\tchrom_pos\tkmer\n", o);
    char ks[32];
    for (uint64_t i = 0; i < total; i++) {
        const KHit *h = &all[i];
        for (int j = 0; j < k; j++)
            ks[j] = CODE_CHARS[(h->code >> (2 * (k - 1 - j))) & 3];
        ks[k] = '\0';
        fprintf(o, "%s\t%u\t%s\t%s\t%" PRIu64 "\t%s\n",
                names[h->seq], h->qpos + 1,
                h->orient ? "RC" : "FWD", g.c[h->chrom].name, h->cpos, ks);
    }
    if (fflush(o) != 0 || ferror(o))
        die("write error on '%s'", pos[3]);
    fclose(o);
    free(all);
    free(slots);
    free(insts);
    for (uint64_t i = 0; i < nnames; i++)
        free(names[i]);
    free(names);
    cgb_close(&g);
}

/* ------------------------------------------------------------------ */
/* help / main                                                         */
/* ------------------------------------------------------------------ */

static void print_help(void)
{
    printf(
"crispr_scan — pack FASTA genomes, scan CRISPR-Cas9 off-target sites, map k-mers\n"
"\n"
"Usage:\n"
"  crispr_scan pack <in.fa> <out.cgb>\n"
"      Pack a plain FASTA genome into the compact .cgb format (CGB1, little-\n"
"      endian: u32 magic 0x43474231, u32 n_chrom, per-chrom records\n"
"      [u32 name_len, name, u64 len], then per-chrom 4-bit-packed blocks,\n"
"      2 bases per byte, low nibble first; A=0 C=1 G=2 T=3, else 4).\n"
"      Case-insensitive; names are the first whitespace-delimited token of\n"
"      the '>' line; blank lines, multi-line sequences and CRLF handled.\n"
"  crispr_scan unpack <bin.cgb> <out.fa>\n"
"      Inverse of pack (code-4 bases are written as 'N'; 60 bases/line).\n"
"  crispr_scan stat <bin.cgb>\n"
"      Print tab-separated lines: \"n_chrom<TAB>N\", \"total_bases<TAB>T\",\n"
"      then one line per chrom in pack order: \"<name><TAB><length>\".\n"
"  crispr_scan scan <genome.cgb> <guides.tsv> <max_mm> <out.tsv> [threads]\n"
"                [--ngg-only] [--canonical-only]\n"
"      Report every genomic site (window + PAM fully inside the chrom)\n"
"      whose window is within max_mm Hamming mismatches of the guide.\n"
"      guides.tsv: first line is a '#' header (skipped), then rows\n"
"      guide_id<TAB>system<TAB>sequence; system is SPCAS9 (20 nt, NGG + NAG\n"
"      PAMs) or SACAS9 (21 nt, NNGRRT + NNGRRV PAMs); A/C/G/T only.\n"
"      --ngg-only restricts SPCAS9 to NGG; --canonical-only restricts\n"
"      SACAS9 to NNGRRT.  Output columns: guide_id, system, chrom, pos\n"
"      (1-based plus-strand window start), strand (+/-), pam (actual PAM\n"
"      bases on the PAM-bearing strand), mismatches, plus_window.  Rows\n"
"      sorted by guide order, mismatches, chrom order, position, strand;\n"
"      byte-identical output for any thread count.\n"
"  crispr_scan kmap <genome.cgb> <query.fa> <k> <out.tsv> [threads]\n"
"      Exact k-mer sharing map, 1 <= k <= 31.  For every ACGT k-mer of\n"
"      every query sequence, report all exact matches on either genome\n"
"      strand.  k-mers with non-ACGT characters are skipped.  Output\n"
"      columns: query_name, query_pos (1-based), orient (FWD|RC), chrom,\n"
"      chrom_pos (1-based), kmer (plus-strand read).\n"
"\n"
"  [threads] defaults to hardware concurrency, capped at 16.\n"
"  Exit status: 0 on success; non-zero with a message on stderr on error.\n");
}

int main(int argc, char **argv)
{
    char_code_init();
    if (argc < 2) {
        fprintf(stderr, "crispr_scan: missing mode (try --help)\n");
        return 1;
    }
    const char *m = argv[1];
    if (strcmp(m, "--help") == 0 || strcmp(m, "-h") == 0) {
        print_help();
        return 0;
    }
    if (strcmp(m, "pack") == 0) {
        if (argc != 4)
            die("pack usage: crispr_scan pack <in.fa> <out.cgb>");
        do_pack(argv[2], argv[3]);
    } else if (strcmp(m, "unpack") == 0) {
        if (argc != 4)
            die("unpack usage: crispr_scan unpack <bin.cgb> <out.fa>");
        do_unpack(argv[2], argv[3]);
    } else if (strcmp(m, "stat") == 0) {
        if (argc != 3)
            die("stat usage: crispr_scan stat <bin.cgb>");
        do_stat(argv[2]);
    } else if (strcmp(m, "scan") == 0) {
        do_scan(argc, argv);
    } else if (strcmp(m, "kmap") == 0) {
        do_kmap(argc, argv);
    } else {
        die("unknown mode '%s' (try --help)", m);
    }
    return 0;
}