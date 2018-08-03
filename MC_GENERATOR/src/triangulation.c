#include <triangulation.h>

/* ************************************************************************************ */
/* ********   Preliminaries - constants, defines, etc  ******************************** */
/* ************************************************************************************ */
/* Various cells (links, triangles, simplices, hypersimplices) are parameterised
   by their vertices in standard hypercube. Complication here is that reflections
   MUST be used to glue neighboring triangulated hypercubes. ... 
   Vertices of each cell are parameterised by one uint number, where groups of
   successive 4 bits represent particular vertex. To ensure that only needed amount
   of bits is present we use the following bit masks, where index stands for the
   cell's dimensionality (number of vertices) */

static const uint cell_mask[6] = {
  0x0,     /* empty */
  0x0F,    /* vertices ( 0-cells ) */
  0xFF,    /* links    ( 1-cells ) */
  0x0FFF,  /* triangles (2-cells with 3 vertices) */
  0xFFFF,  /* simplices (3-cells with 4 vertices) */
  0x0FFFFF /* hypersimplices (4-sim simplices with 5 vertices ) */
};
/* ------------------------------------------------------------------------- */
/* Vertices of hypersimplices for parity = 0. For non canonical parity use reflection
   (see below) to get bit mask of hypersimplices */
static const uint hs_v_table[HYPERSIMPLEX_PER_SITE] = {
  0x95310, 0xa6320, 0xc6540, 0xca980, 0xfba93, 0xfdc95, 0xfeca6, 0xf7653,
  0xf9530, 0xfa630, 0xfc650, 0xfca90, 0xfa930, 0xfc950, 0xfca60, 0xf6530
};
/* ------------------------------------------------------------------------- */
/* Links enumeration. Link with mask v_mask is 'my' iff
    (v_mask has no 'common bits')
    || (all 'common bits' in v_mask [maximum 3] are zero)
   Enumeration of 'my' links: [parity][k] = {vertex mask}
   Masks below are ascending to speed up searches. */
static const uint l_enum_table[MAX_PARITY][LINK_PER_SITE] = {
  { 0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x80, 0x90, 0xa0, 0xc0, 0xf0 },
  { 0x01, 0x20, 0x21, 0x40, 0x41, 0x42, 0x80, 0x81, 0x82, 0x84, 0xe1 },
  { 0x02, 0x10, 0x12, 0x40, 0x41, 0x42, 0x80, 0x81, 0x82, 0x84, 0xd2 },
  { 0x01, 0x02, 0x03, 0x40, 0x50, 0x60, 0x80, 0x90, 0xa0, 0xc0, 0xc3 },
  { 0x04, 0x10, 0x14, 0x20, 0x21, 0x24, 0x80, 0x81, 0x82, 0x84, 0xb4 },
  { 0x01, 0x04, 0x05, 0x06, 0x20, 0x30, 0x80, 0x90, 0xa0, 0xa5, 0xc0 },
  { 0x02, 0x03, 0x04, 0x05, 0x06, 0x10, 0x80, 0x90, 0x96, 0xa0, 0xc0 },
  { 0x01, 0x02, 0x04, 0x12, 0x14, 0x24, 0x80, 0x81, 0x82, 0x84, 0x87 },
  { 0x08, 0x10, 0x18, 0x20, 0x21, 0x28, 0x40, 0x41, 0x42, 0x48, 0x78 },
  { 0x01, 0x08, 0x09, 0x0a, 0x0c, 0x20, 0x30, 0x40, 0x50, 0x60, 0x69 },
  { 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0c, 0x10, 0x40, 0x50, 0x5a, 0x60 },
  { 0x01, 0x02, 0x08, 0x12, 0x18, 0x28, 0x40, 0x41, 0x42, 0x48, 0x4b },
  { 0x04, 0x05, 0x06, 0x08, 0x09, 0x0a, 0x0c, 0x10, 0x20, 0x30, 0x3c },
  { 0x01, 0x04, 0x08, 0x14, 0x18, 0x20, 0x21, 0x24, 0x28, 0x2d, 0x48 },
  { 0x02, 0x04, 0x08, 0x10, 0x12, 0x14, 0x18, 0x1e, 0x24, 0x28, 0x48 },
  { 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x08, 0x09, 0x0a, 0x0c, 0x0f }
};
/* ------------------------------------------------------------------------- */
/* Triangles enumeration. Triangle with mask v_mask is 'my' iff
    (v_mask has no 'common bits')
    || (all 'common bits' in v_mask [maximum 2] are zero)
   Enumeration of 'my' triangles:  [parity][k] = {vertex mask}
   Masks below are ascending to speed up searches. */
static const uint tr_enum_table[MAX_PARITY][TRIANGLE_PER_SITE] = {
  { 0x310, 0x320, 0x510, 0x530, 0x540, 0x620, 0x630, 0x640, 0x650, 0x653, 0x910, 0x930, 0x950, 0x980, 0xa20, 0xa30, 0xa60,
    0xa80, 0xa90, 0xa93, 0xc40, 0xc50, 0xc60, 0xc80, 0xc90, 0xc95, 0xca0, 0xca6, 0xf30, 0xf50, 0xf60, 0xf90, 0xfa0, 0xfc0 },
  { 0x201, 0x231, 0x401, 0x420, 0x421, 0x451, 0x642, 0x721, 0x741, 0x742, 0x801, 0x820, 0x821, 0x840, 0x841, 0x842, 0x891,
    0xa82, 0xb21, 0xb81, 0xb82, 0xc84, 0xd41, 0xd81, 0xd84, 0xe21, 0xe41, 0xe42, 0xe71, 0xe81, 0xe82, 0xe84, 0xeb1, 0xed1 },
  { 0x102, 0x132, 0x402, 0x410, 0x412, 0x462, 0x471, 0x472, 0x541, 0x712, 0x802, 0x810, 0x812, 0x840, 0x841, 0x842, 0x8a2,
    0x8b1, 0x8b2, 0x981, 0xb12, 0xc84, 0xd12, 0xd41, 0xd42, 0xd72, 0xd81, 0xd82, 0xd84, 0xdb2, 0xde2, 0xe42, 0xe82, 0xe84 },
  { 0x013, 0x023, 0x450, 0x460, 0x501, 0x503, 0x560, 0x563, 0x602, 0x603, 0x890, 0x8a0, 0x901, 0x903, 0x950, 0x9a0, 0x9a3,
    0xa02, 0xa03, 0xa60, 0xc03, 0xc40, 0xc50, 0xc53, 0xc60, 0xc63, 0xc80, 0xc90, 0xc93, 0xc95, 0xca0, 0xca3, 0xca6, 0xcf3 },
  { 0x104, 0x154, 0x174, 0x204, 0x210, 0x214, 0x217, 0x264, 0x274, 0x321, 0x804, 0x810, 0x814, 0x820, 0x821, 0x824, 0x8c4,
    0x8d1, 0x8d4, 0x8e2, 0x8e4, 0x981, 0xa82, 0xb14, 0xb21, 0xb24, 0xb74, 0xb81, 0xb82, 0xb84, 0xbd4, 0xbe4, 0xd14, 0xe24 },
  { 0x015, 0x045, 0x064, 0x065, 0x206, 0x230, 0x301, 0x305, 0x306, 0x365, 0x890, 0x8c0, 0x901, 0x905, 0x930, 0x9c0, 0x9c5,
    0xa05, 0xa06, 0xa20, 0xa30, 0xa35, 0xa65, 0xa80, 0xa90, 0xa93, 0xa95, 0xac0, 0xac5, 0xac6, 0xaf5, 0xc04, 0xc05, 0xc06 },
  { 0x026, 0x032, 0x035, 0x036, 0x046, 0x054, 0x056, 0x103, 0x105, 0x356, 0x8a0, 0x8c0, 0x903, 0x905, 0x906, 0x910, 0x936,
    0x956, 0x980, 0x9a0, 0x9a3, 0x9a6, 0x9c0, 0x9c5, 0x9c6, 0x9f6, 0xa02, 0xa03, 0xa06, 0xac0, 0xac6, 0xc04, 0xc05, 0xc06 },
  { 0x012, 0x014, 0x024, 0x123, 0x124, 0x127, 0x145, 0x147, 0x246, 0x247, 0x801, 0x802, 0x804, 0x812, 0x814, 0x817, 0x824, 
    0x827, 0x847, 0x891, 0x8a2, 0x8b1, 0x8b2, 0x8b7, 0x8c4, 0x8d1, 0x8d4, 0x8d7, 0x8e2, 0x8e4, 0x8e7, 0xb12, 0xd14, 0xe24 },
  { 0x108, 0x198, 0x1b8, 0x1d8, 0x208, 0x210, 0x218, 0x21b, 0x2a8, 0x2b8, 0x2e8, 0x321, 0x408, 0x410, 0x418, 0x41d, 0x420,
    0x421, 0x428, 0x42e, 0x4c8, 0x4d8, 0x4e8, 0x541, 0x642, 0x718, 0x721, 0x728, 0x741, 0x742, 0x748, 0x7b8, 0x7d8, 0x7e8 },
  { 0x019, 0x089, 0x0a8, 0x0a9, 0x0c8, 0x0c9, 0x0ca, 0x20a, 0x230, 0x301, 0x309, 0x30a, 0x3a9, 0x40c, 0x450, 0x501, 0x509,
    0x50c, 0x530, 0x5c9, 0x609, 0x60a, 0x60c, 0x620, 0x630, 0x639, 0x640, 0x650, 0x653, 0x659, 0x6a9, 0x6c9, 0x6ca, 0x6f9 },
  { 0x02a, 0x032, 0x039, 0x03a, 0x08a, 0x098, 0x09a, 0x0c8, 0x0c9, 0x0ca, 0x103, 0x109, 0x39a, 0x40c, 0x460, 0x503, 0x509,
    0x50a, 0x50c, 0x510, 0x53a, 0x540, 0x560, 0x563, 0x56a, 0x59a, 0x5c9, 0x5ca, 0x5fa, 0x602, 0x603, 0x60a, 0x60c, 0x6ca },
  { 0x012, 0x018, 0x028, 0x123, 0x128, 0x12b, 0x189, 0x18b, 0x1d8, 0x28a, 0x28b, 0x2e8, 0x401, 0x402, 0x408, 0x412, 0x418,
    0x41b, 0x41d, 0x428, 0x42b, 0x42e, 0x451, 0x462, 0x471, 0x472, 0x47b, 0x48b, 0x4c8, 0x4d8, 0x4db, 0x4e8, 0x4eb, 0x712 },
  { 0x04c, 0x054, 0x059, 0x05c, 0x064, 0x065, 0x06a, 0x06c, 0x08c, 0x098, 0x09c, 0x0a8, 0x0a9, 0x0ac, 0x105, 0x109, 0x206,
    0x20a, 0x305, 0x306, 0x309, 0x30a, 0x30c, 0x310, 0x320, 0x35c, 0x365, 0x36c, 0x39c, 0x3a9, 0x3ac, 0x3fc, 0x59c, 0x6ac },
  { 0x014, 0x018, 0x048, 0x145, 0x148, 0x14d, 0x174, 0x189, 0x18d, 0x1b8, 0x201, 0x204, 0x208, 0x214, 0x217, 0x218, 0x21b,
    0x21d, 0x231, 0x248, 0x24d, 0x24e, 0x264, 0x274, 0x27d, 0x28d, 0x28e, 0x2a8, 0x2b8, 0x2bd, 0x2ed, 0x48c, 0x48d, 0x48e },
  { 0x024, 0x028, 0x048, 0x102, 0x104, 0x108, 0x124, 0x127, 0x128, 0x12b, 0x12e, 0x132, 0x147, 0x148, 0x14d, 0x14e, 0x154,
    0x17e, 0x18b, 0x18d, 0x18e, 0x198, 0x1be, 0x1de, 0x246, 0x247, 0x248, 0x24e, 0x28a, 0x28b, 0x28e, 0x48c, 0x48d, 0x48e },
  { 0x013, 0x015, 0x019, 0x023, 0x026, 0x02a, 0x035, 0x036, 0x039, 0x03a, 0x03f, 0x045, 0x046, 0x04c, 0x056, 0x059, 0x05c,
    0x05f, 0x06a, 0x06c, 0x06f, 0x089, 0x08a, 0x08c, 0x09a, 0x09c, 0x09f, 0x0ac, 0x0af, 0x0cf, 0x356, 0x39a, 0x59c, 0x6ac }
};
/* ------------------------------------------------------------------------- */
/* Simplex enumeration. Simplex with mask v_mask is 'my' iff
    (v_mask has no 'common bits')
    || ('common bit' in v_mask [only one is possible] is zero)
   Enumeration of 'my' simplices:  [parity][k] = {vertex mask}
   Masks below are ascending to speed up searches. */
static const uint s_enum_table[MAX_PARITY][SIMPLEX_PER_SITE] = {
  { 0x5310, 0x6320, 0x6530, 0x6540, 0x7653, 0x9310, 0x9510, 0x9530, 0xa320, 0xa620,
    0xa630, 0xa930, 0xa980, 0xba93, 0xc540, 0xc640, 0xc650, 0xc950, 0xc980, 0xca60,
    0xca80, 0xca90, 0xdc95, 0xeca6, 0xf530, 0xf630, 0xf650, 0xf653, 0xf930, 0xf950,
    0xfa30, 0xfa60, 0xfa90, 0xfa93, 0xfc50, 0xfc60, 0xfc90, 0xfc95, 0xfca0, 0xfca6 },
  { 0x4201, 0x6742, 0x7231, 0x7421, 0x7451, 0x8201, 0x8401, 0x8420, 0x8421, 0xab82,
    0xb231, 0xb721, 0xb821, 0xb891, 0xcd84, 0xd451, 0xd741, 0xd841, 0xd891, 0xdb81,
    0xe421, 0xe642, 0xe721, 0xe741, 0xe742, 0xe821, 0xe841, 0xe842, 0xea82, 0xeb21,
    0xeb71, 0xeb81, 0xeb82, 0xec84, 0xed41, 0xed71, 0xed81, 0xed84, 0xedb1, 0xedb7 },
  { 0x4102, 0x4712, 0x4762, 0x5471, 0x7132, 0x8102, 0x8402, 0x8410, 0x8412, 0x8b12,
    0x8ba2, 0x98b1, 0xb132, 0xb712, 0xce84, 0xd412, 0xd471, 0xd472, 0xd541, 0xd712,
    0xd812, 0xd841, 0xd842, 0xd8b1, 0xd8b2, 0xd981, 0xdb12, 0xdb72, 0xdc84, 0xde42,
    0xde72, 0xde82, 0xde84, 0xdeb2, 0xdeb7, 0xe462, 0xe472, 0xe842, 0xe8a2, 0xe8b2 },
  { 0x4560, 0x5013, 0x5603, 0x5673, 0x6023, 0x89a0, 0x9013, 0x9501, 0x9503, 0x9a03,
    0x9ab3, 0xa023, 0xa602, 0xa603, 0xc450, 0xc460, 0xc503, 0xc560, 0xc563, 0xc603,
    0xc890, 0xc8a0, 0xc903, 0xc950, 0xc953, 0xc9a0, 0xc9a3, 0xca03, 0xca60, 0xca63,
    0xcd95, 0xcea6, 0xcf53, 0xcf63, 0xcf93, 0xcf95, 0xcfa3, 0xcfa6, 0xf563, 0xf9a3 },
  { 0x1754, 0x2104, 0x2174, 0x2764, 0x3217, 0x8104, 0x8204, 0x8210, 0x8214, 0x8d14,
    0x8dc4, 0x8e24, 0x8ec4, 0x8ed4, 0x98d1, 0xa8e2, 0xb174, 0xb214, 0xb217, 0xb274,
    0xb321, 0xb814, 0xb821, 0xb824, 0xb8d1, 0xb8d4, 0xb8e2, 0xb8e4, 0xb981, 0xba82,
    0xbd14, 0xbd74, 0xbe24, 0xbe74, 0xbed4, 0xbed7, 0xd154, 0xd174, 0xe264, 0xe274 },
  { 0x0645, 0x2306, 0x3015, 0x3065, 0x3675, 0x89c0, 0x9015, 0x9301, 0x9305, 0x9c05,
    0x9cd5, 0x9fc5, 0xa065, 0xa206, 0xa230, 0xa305, 0xa306, 0xa365, 0xa890, 0xa8c0,
    0xa905, 0xa930, 0xa935, 0xa9c0, 0xa9c5, 0xa9f3, 0xa9f5, 0xab93, 0xac05, 0xac06,
    0xac65, 0xaec6, 0xaf35, 0xaf65, 0xafc5, 0xafc6, 0xc045, 0xc064, 0xc065, 0xf365 },
  { 0x0326, 0x0356, 0x0546, 0x1035, 0x3576, 0x8ac0, 0x9035, 0x9036, 0x9056, 0x9103,
    0x9105, 0x9356, 0x98a0, 0x98c0, 0x9a03, 0x9a06, 0x9a36, 0x9ac0, 0x9ac6, 0x9af3,
    0x9af6, 0x9ba3, 0x9c05, 0x9c06, 0x9c56, 0x9cf5, 0x9cf6, 0x9dc5, 0x9f36, 0x9f56,
    0xa026, 0xa032, 0xa036, 0xac06, 0xace6, 0xacf6, 0xc046, 0xc054, 0xc056, 0xf356 },
  { 0x0124, 0x1237, 0x1247, 0x1457, 0x2467, 0x8012, 0x8014, 0x8024, 0x8124, 0x8127,
    0x8147, 0x8247, 0x89b1, 0x89d1, 0x8ab2, 0x8ae2, 0x8b12, 0x8b17, 0x8b27, 0x8bd1,
    0x8bd7, 0x8be2, 0x8be7, 0x8cd4, 0x8ce4, 0x8d14, 0x8d17, 0x8d47, 0x8de4, 0x8de7,
    0x8e24, 0x8e27, 0x8e47, 0xb123, 0xb127, 0xbde7, 0xd145, 0xd147, 0xe246, 0xe247 },
  { 0x1b98, 0x1d98, 0x1db8, 0x2108, 0x21b8, 0x2ba8, 0x2ea8, 0x2eb8, 0x321b, 0x4108,
    0x41d8, 0x4208, 0x4210, 0x4218, 0x42e8, 0x4dc8, 0x4ec8, 0x4ed8, 0x541d, 0x642e,
    0x71b8, 0x71d8, 0x7218, 0x721b, 0x72b8, 0x72e8, 0x7321, 0x7418, 0x741d, 0x7421,
    0x7428, 0x742e, 0x74d8, 0x74e8, 0x7541, 0x7642, 0x7db8, 0x7eb8, 0x7ed8, 0x7edb },
  { 0x0a89, 0x0c89, 0x0ca8, 0x0ca9, 0x230a, 0x3019, 0x30a9, 0x3ab9, 0x3fa9, 0x450c,
    0x5019, 0x50c9, 0x5301, 0x5309, 0x5cd9, 0x5fc9, 0x60a9, 0x60c9, 0x60ca, 0x620a,
    0x6230, 0x6309, 0x630a, 0x63a9, 0x63f9, 0x640c, 0x6450, 0x6509, 0x650c, 0x6530,
    0x6539, 0x653f, 0x65c9, 0x65f9, 0x6753, 0x6ca9, 0x6eca, 0x6fa9, 0x6fc9, 0x6fca },
  { 0x032a, 0x039a, 0x098a, 0x0c8a, 0x0c98, 0x0c9a, 0x1039, 0x39ba, 0x3f9a, 0x460c,
    0x5039, 0x503a, 0x509a, 0x50c9, 0x50ca, 0x5103, 0x5109, 0x539a, 0x53fa, 0x540c,
    0x5460, 0x5603, 0x560a, 0x560c, 0x563a, 0x563f, 0x56ca, 0x56fa, 0x5763, 0x5c9a,
    0x5cf9, 0x5cfa, 0x5dc9, 0x5f9a, 0x602a, 0x6032, 0x603a, 0x60ca, 0x6cea, 0x6cfa },
  { 0x0128, 0x123b, 0x128b, 0x189b, 0x1d89, 0x1d8b, 0x28ab, 0x2e8a, 0x2e8b, 0x4012,
    0x4018, 0x4028, 0x4128, 0x412b, 0x418b, 0x41d8, 0x41db, 0x428b, 0x42e8, 0x42eb,
    0x451d, 0x4571, 0x462e, 0x4672, 0x4712, 0x471b, 0x471d, 0x472b, 0x472e, 0x47db,
    0x47eb, 0x4cd8, 0x4ce8, 0x4d8b, 0x4de8, 0x4deb, 0x4e8b, 0x7123, 0x712b, 0x7deb },
  { 0x054c, 0x059c, 0x064c, 0x0654, 0x065c, 0x06ac, 0x098c, 0x0a8c, 0x0a98, 0x0a9c,
    0x1059, 0x206a, 0x3059, 0x305c, 0x3065, 0x306a, 0x306c, 0x309c, 0x30a9, 0x30ac,
    0x3105, 0x3109, 0x3206, 0x320a, 0x359c, 0x35fc, 0x365c, 0x365f, 0x36ac, 0x36fc,
    0x3765, 0x39fc, 0x3a9c, 0x3a9f, 0x3afc, 0x3ba9, 0x59dc, 0x59fc, 0x6aec, 0x6afc },
  { 0x0148, 0x145d, 0x148d, 0x1745, 0x174d, 0x189d, 0x1b89, 0x1b8d, 0x2014, 0x2018,
    0x2048, 0x2148, 0x214d, 0x2174, 0x217b, 0x217d, 0x218d, 0x21b8, 0x21bd, 0x2317,
    0x231b, 0x248d, 0x248e, 0x24ed, 0x264e, 0x2674, 0x274d, 0x274e, 0x27bd, 0x27ed,
    0x28ed, 0x2a8e, 0x2ab8, 0x2b8d, 0x2b8e, 0x2bed, 0x48cd, 0x48ec, 0x48ed, 0x7bed },
  { 0x0248, 0x1024, 0x1028, 0x1048, 0x1247, 0x1248, 0x124e, 0x127b, 0x127e, 0x128b,
    0x128e, 0x12be, 0x1327, 0x132b, 0x147d, 0x147e, 0x148d, 0x148e, 0x14de, 0x1547,
    0x154d, 0x17be, 0x17de, 0x18bd, 0x18be, 0x18de, 0x198b, 0x198d, 0x1bde, 0x246e,
    0x2476, 0x247e, 0x248e, 0x28ae, 0x28ba, 0x28be, 0x48ce, 0x48dc, 0x48de, 0x7bde },
  { 0x0135, 0x0139, 0x0159, 0x0236, 0x023a, 0x026a, 0x0356, 0x0359, 0x035f, 0x036a,
    0x036f, 0x039a, 0x039f, 0x03af, 0x0456, 0x045c, 0x046c, 0x056c, 0x056f, 0x059c,
    0x059f, 0x05cf, 0x06ac, 0x06af, 0x06cf, 0x089a, 0x089c, 0x08ac, 0x09ac, 0x09af,
    0x09cf, 0x0acf, 0x3567, 0x356f, 0x39ab, 0x39af, 0x59cd, 0x59cf, 0x6ace, 0x6acf }
};
/* ************************************************************************************ */
/* ******    Cell dependencies   ****************************************************** */
/* ************************************************************************************ */
/* How many different cells depend on a particular site (evidently, this number varies with
   site parity). Various cells are represented by their dimensionality. */
static const uint site_dependence_dim[4][MAX_PARITY] = {
  { 48, 8, 8, 32, 8, 32, 32, 8, 8, 32, 32, 8, 32, 8, 8, 48 }, /* links */
  { 240, 24, 24, 160, 24, 160, 160, 24, 24, 160, 160, 24, 160, 24, 24, 240 }, /* triangles */
  { 384, 32, 32, 256, 32, 256, 256, 32, 32, 256, 256, 32, 256, 32, 32, 384 }, /* simplices */
  { 192, 16, 16, 128, 16, 128, 128, 16, 16, 128, 128, 16, 128, 16, 16, 192 } /* hypersimplices */
};

/*
  Road map from sites to cells of various dimension, constructed dynamically.
  Entries are 0xZSII (hex), where:
  II (8 bits) - local index of target cell;
  S  (4 bits) - shift mask, how to get from site to target cell;
  Z  (4 bits) - shift sign, if bit is set than move downwards;
 */
static uint *site_dependence[4][MAX_PARITY];
static uchar site_dependence_done = 0;

/* ---------------------------------------------------------------------------------------- */
/* How many different triangles depend on a particular link (this number varies with parity
   and 'small' link index).*/
static const uchar link2triangle_dim[MAX_PARITY][LINK_PER_SITE] = {
  { 6, 6, 14, 6, 14, 14, 6, 14, 14, 14, 6 }, { 6, 6, 14, 6, 14, 10, 6, 14, 10, 10, 6 },
  { 6, 6, 14, 6, 10, 14, 6, 10, 14, 10, 6 }, { 6, 6, 14, 6, 10, 10, 6, 10, 10, 14, 6 },
  { 6, 6, 14, 6, 10, 14, 6, 10, 10, 14, 6 }, { 6, 6, 14, 10, 6, 10, 6, 10, 14, 6, 10 },
  { 6, 10, 6, 10, 14, 6, 6, 14, 6, 10, 10 }, { 6, 6, 6, 10, 10, 10, 6, 14, 14, 14, 6 },
  { 6, 6, 14, 6, 10, 14, 6, 10, 10, 14, 6 }, { 6, 6, 14, 10, 10, 6, 10, 6, 10, 14, 6 },
  { 6, 10, 6, 10, 14, 10, 6, 6, 14, 6, 10 }, { 6, 6, 6, 10, 10, 10, 6, 14, 14, 14, 6 },
  { 6, 10, 10, 6, 10, 10, 14, 6, 6, 14, 6 }, { 6, 6, 6, 10, 10, 6, 14, 14, 14, 6, 10 },
  { 6, 6, 6, 6, 14, 14, 14, 6, 10, 10, 10 }, { 6, 6, 14, 6, 14, 14, 6, 14, 14, 14, 6 }
};
/*
  Road map from links to triangles constructed dynamically. Entries are 0xGZSII (hex), where:
  II (8 bits) - local index of target triangle;
  S  (4 bits) - shift mask, how to get from link 'base point' to target 'base';
  Z  (4 bits) - shift sign, if bit is set than move downwards;
  G  (4 bits) - if non-zero then link enters with negative sign to triangle.
 */
static uint *link2triangle[MAX_PARITY][LINK_PER_SITE];
static uchar link2triangle_done = 0;


#define  INDEX_MASK       0x000FF
#define  SHIFT_MASK       0x00F00
#define  SHIFT_SIGN       0x0F000
#define  SIGN_MASK        0xF0000
/* ************************************************************************************ */
/* *******  Helping routines  ********************************************************* */
/* ************************************************************************************ */
/* Reflect bits in mask according to parity */
static uint reflection( uint mask, uchar parity, uchar dim ){
  uchar i;
  parity |= (parity & 0x0F) << 4;
  for( i = 0 ; i < 4; i++ ) mask ^= parity << (8*i);
  return mask & cell_mask[dim];
};
/* ------------------------------------------------------------------------- */
/* Erases particular vertex pointed by 'where' from vertex mask 'mask' */
static uint erase_bits( uint mask, uchar where , uchar dim ){
  uchar k;
  uint low = 0, high = 0;
  for( k = 0; k < where; k++ ) low |= 0x0F << (4*k);
  for( k = where+1; k < dim; k++ ) high |= 0x0F << (4*k);
  return (((mask & high) >> 4) | (mask & low)) & cell_mask[dim-1];
};
/* ------------------------------------------------------------------------- */
/* Shifts reference index m according to bit mask, in which only 4 lowest bits
   are significant. Shifts upwards/downwards depending on sign. */
static uint shift_index( uint m, uchar v_mask, int sign ){
  uchar i;
  for( i = 0; i < 4; i++ )
    if( v_mask & (1<<i) )
      m = ( sign > 0 ) ? index_up( m, i ) : index_down( m, i );
  return m;
};
/* ------------------------------------------------------------------------- */
static char get_parity( uint m , uchar *parity ){
  uchar i, x[4];
  char ret = 1;
  site_coordinates( x, m );
  for( *parity = 0, i = 0; i < 4; i++ )
    if( x[i] % 2 ){ *parity |= 1 << i; ret *= -1; };
  return ret;
};
/* ------------------------------------------------------------------------- */
/* Searches and counts common bits present in bit mask 'num'. 'num' represents
   vertex mask of links, triangles, simplices (which is accounted for by 'dim'
   argument). E.g. link mask 0x10 has four common bits, all of which are zero.
   Output:  where - location of common bits in vertex mask; which - which bits
   are common. */
static uchar common_bits( uint num, uchar *where, uchar *which, uchar dim ){
  uint mask = 0x11111111 & cell_mask[dim];
  uchar i, ret = 0;
  for( i = 0; i < 4; i++, mask <<= 1 ){
    if( mask == (num & mask) ){
      if( where ) where[ret] = i;
      if( which ) which[ret] = 1;
      ret++;
    };
    if( mask == (~num & mask) ){
      if( where ) where[ret] = i;
      if( which ) which[ret] = 0;
      ret++;
    };
  };
#ifdef DEBUG
  if( ret > (5 - dim) ) panic("Too many common bits");
#endif
  return ret;
};
/* ------------------------------------------------------------------------- */
static inline uchar l_by_mask( uint mask, uchar parity, uchar start, uchar end ){
  uchar middle = (start + end)/2;
  if( mask == l_enum_table[parity][start] ) return start;
#ifdef DEBUG
  if( start >= end ) panic2("Link mask 0x%x was not found", mask );
#endif
  if( mask == l_enum_table[parity][middle] ) return middle;
  if( mask >  l_enum_table[parity][middle] ) return l_by_mask( mask, parity, middle, end );
  return l_by_mask( mask, parity, start, middle);
};
static uint l_enum( uint m, uchar parity, uint mask ){
  uchar i, where[3], which[3];
  uchar cb = common_bits( mask, where, which, 2 );
  if( (cb == 1 && which[0])
      || ( cb == 2 && (which[0] || which[1]) ) 
      || ( cb == 3 && (which[0] || which[1] || which[2]) ) ){ /* NOT my link */
    for( i = 0; i < cb ; i++ ){
      if( !which[i] ) continue;
      parity ^= 1 << where[i];
      mask = reflection( mask, 1 << where[i], 2 );
      m = index_up( m, where[i] );
    };
    return l_enum( m, parity, mask );
  };
  i = l_by_mask( mask, parity, 0, LINK_PER_SITE );
  return L_ENUM( i , m );
};
/* ------------------------------------------------------------------------- */
static inline uchar tr_by_mask( uint mask, uchar parity, uchar start, uchar end ){
  uchar middle = (start + end)/2;
  if( mask == tr_enum_table[parity][start] ) return start;
#ifdef DEBUG
  if( start >= end ) panic2("Triangle mask 0x%x was not found", mask );
#endif
  if( mask == tr_enum_table[parity][middle] ) return middle;
  if( mask > tr_enum_table[parity][middle] ) return tr_by_mask( mask, parity, middle, end );
  return tr_by_mask( mask, parity, start, middle);
};
static uint tr_enum( uint m, uchar parity, uint mask ){
  uchar i, where[2], which[2];
  uchar cb = common_bits( mask, where, which, 3 );
  if( (cb == 1 && which[0])
      || ( cb == 2 && (which[0] || which[1]) ) ){ /* NOT my triangle */
    for( i = 0; i < cb ; i++ ){
      if( !which[i] ) continue;
      parity ^= 1 << where[i];
      mask = reflection( mask, 1 << where[i], 3 );
      m = index_up( m, where[i] );
    };
    return tr_enum( m, parity, mask );
  };
  i = tr_by_mask( mask, parity, 0, TRIANGLE_PER_SITE );
  return TR_ENUM( i , m );
};
/* ------------------------------------------------------------------------- */
static inline uchar s_by_mask( uint mask, uchar parity, uchar start, uchar end ){
  uchar middle = (start + end)/2;
  if( mask == s_enum_table[parity][start] ) return start;
#ifdef DEBUG
  if( start >= end ) panic2("Simplex mask 0x%x was not found", mask );
#endif
  if( mask == s_enum_table[parity][middle] ) return middle;
  if( mask > s_enum_table[parity][middle] ) return s_by_mask( mask, parity, middle, end );
  return s_by_mask( mask, parity, start, middle);
};
static uint s_enum( uint m, uchar parity, uint mask ){
  uchar i, where, which;
  uchar cb = common_bits( mask, &where, &which, 4 );
  if( cb == 1 && which ){ /* NOT my simplex */
    parity ^= 1 << where;
    mask = reflection( mask, 1 << where, 4 );
    m = index_up( m, where );
    return s_enum( m, parity, mask );
  };
  i = s_by_mask( mask, parity, 0, SIMPLEX_PER_SITE );
  return S_ENUM( i , m );
};
/* ************************************************************************************ */
/* ******   Differentials from/to various dimensions   ******************************** */
/* ************************************************************************************ */
/* Differential from simplices to hypersimplices.
   Input:  uint hs - index of hypersimplex on which the differential is to be evaluated;
   Output: uint *s_idx, int *s_sign - five indices and signs of simplices which build up
   the differential:    (d)_hs = \sum sign_k simplex_k */
void Tr_diff_simplex( uint hs, uint *s_idx, int *s_sign ){
  uint m, i_hs;
  uchar i, parity ;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( hs >= HYPERSIMPLEX_PER_SITE * param->ipw[param->D] ) panic("Illegal parameter");
  m = hs/((uint)HYPERSIMPLEX_PER_SITE);
  i_hs = hs - HYPERSIMPLEX_PER_SITE * m;
  get_parity( m, &parity );
  for( i = 0; i < 5; i++ ){
    s_idx[i] = s_enum( m, parity, reflection( erase_bits( hs_v_table[i_hs], i, 5 ), parity, 4 ) );
    s_sign[i] = (i % 2) ? -1 : 1;
  };
};
/* ------------------------------------------------------------------------- */
/* Differential from triangles to simplices.
   Input:  uint s - index of simplex on which the differential is to be evaluated;
   Output: uint *tr_idx, int *tr_sign - four indices and signs of triangles which build up
   the differential:    (d)_s = \sum sign_k triangle_k */
void Tr_diff_triangle( uint s, uint *tr_idx, int *tr_sign ){
  uint m, i_s;
  uchar i, parity ;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( s >= SIMPLEX_PER_SITE * param->ipw[param->D] ) panic("Illegal parameter");
  m = s/((uint)SIMPLEX_PER_SITE);
  i_s = s - SIMPLEX_PER_SITE * m;
  get_parity( m, &parity );
  for( i = 0; i < 4; i++ ){
    tr_idx[i] = tr_enum( m, parity, erase_bits( s_enum_table[parity][i_s], i, 4 ) );
    tr_sign[i] = (i % 2) ? -1 : 1;
  };
};
/* ------------------------------------------------------------------------- */
/* Differential from links to triangles.
   Input:  uint tr - index of triangle on which the differential is to be evaluated;
   Output: uint *l_idx, int *l_sign - three indices and signs of links which build up
   the differential:    (d)_tr = \sum sign_k link_k */
void Tr_diff_link( uint tr, uint *l_idx, int *l_sign ){
  uint m, i_tr;
  uchar i, parity ;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( tr >= TRIANGLE_PER_SITE * param->ipw[param->D] ) panic("Illegal parameter");
  m = tr/((uint)TRIANGLE_PER_SITE);
  i_tr = tr - TRIANGLE_PER_SITE * m;
  get_parity( m, &parity );
  for( i = 0; i < 3; i++ ){
    l_idx[i] = l_enum( m, parity, erase_bits( tr_enum_table[parity][i_tr], i, 3 )  );
    l_sign[i] = (i % 2) ? -1 : 1;
  };
};
/* ------------------------------------------------------------------------- */
/* Differential from sites to links.
   Input:  uint l - index of link on which the differential is to be evaluated;
   Output: uint *s_idx, int *s_sign - two indices and signs of sites which build up
   the differential:    (d)_l = \sum sign_k site_k */
void Tr_diff_site( uint l, uint *idx, int *sign ){
  uint m, i_l;
  uchar i, parity ;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( l >= LINK_PER_SITE * param->ipw[param->D] ) panic("Illegal parameter");
  m = l/((uint)LINK_PER_SITE);
  i_l = l - LINK_PER_SITE * m;
  get_parity( m, &parity );
  for( i = 0; i < 2; i++ ){
    idx[i] = shift_index( m, erase_bits( l_enum_table[parity][i_l], i, 2 ), 1 );
    sign[i] = (i % 2) ? -1 : 1;
  };
};
/* ------------------------------------------------------------------------- */
/* Check the nilpotency of the above differentials */
static void Tr_check_diff(){
  uint m, n, idx[5];
  int sign[5];
  double *st = NULL, *lnk = NULL, *tr = NULL, *s = NULL, *hs = NULL;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( m = 0; m < 4; m++ )
    if( param->size[m] % 2 ) panic("Lattice MUST be even");
  if( !(st = (double *) calloc( param->ipw[param->D], sizeof(double)))
      || !(lnk = (double *) calloc( LINK_PER_SITE * param->ipw[param->D], sizeof(double)))
      || !(tr = (double *) calloc( TRIANGLE_PER_SITE * param->ipw[param->D], sizeof(double)))
      || !(s = (double *) calloc( SIMPLEX_PER_SITE * param->ipw[param->D], sizeof(double)))
      || !(hs = (double *) calloc( HYPERSIMPLEX_PER_SITE * param->ipw[param->D], sizeof(double))) )
    panic("Memory");
  /* ======= diff sites --> links --> triangles =========== */
  for( m = 0; m < param->ipw[param->D]; m++ ) st[m] = 2.0 * RND() - 1.0;
  for( m = 0; m < LINK_PER_SITE * param->ipw[param->D]; m++ ){
    Tr_diff_site( m, idx, sign );
    for( lnk[m] = 0.0, n = 0; n < 2; n++ ) lnk[m] += sign[n] * st[idx[n]];
  };
  for( m = 0; m < TRIANGLE_PER_SITE * param->ipw[param->D]; m++ ){
    Tr_diff_link( m, idx, sign );
    for( tr[m] = 0.0, n = 0; n < 3; n++ ) tr[m] += sign[n] * lnk[idx[n]];
    if( fabs(tr[m]) > FLT_EPSILON )
      panic("Diff sites --> links --> triangles FAILED");
  };
  /* ======= diff links --> triangles --> simplices =========== */
  for( m = 0; m < LINK_PER_SITE * param->ipw[param->D]; m++ ) lnk[m] = 2.0 * RND() - 1.0;
  for( m = 0; m < TRIANGLE_PER_SITE * param->ipw[param->D]; m++ ){
    Tr_diff_link( m, idx, sign );
    for( tr[m] = 0.0, n = 0; n < 3; n++ ) tr[m] += sign[n] * lnk[idx[n]];
  };
  for( m = 0; m < SIMPLEX_PER_SITE * param->ipw[param->D]; m++ ){
    Tr_diff_triangle( m, idx, sign );
    for( s[m] = 0.0, n = 0; n < 4; n++ ) s[m] += sign[n] * tr[idx[n]];
    if( fabs(s[m]) > FLT_EPSILON )
      panic("Diff links --> triangles --> simplex FAILED");
  };
  /* ======= diff triangles --> simplices --> hypersimplices =========== */
  for( m = 0; m < TRIANGLE_PER_SITE * param->ipw[param->D]; m++ ) tr[m] = 2.0 * RND() - 1.0;
  for( m = 0; m < SIMPLEX_PER_SITE * param->ipw[param->D]; m++ ){
    Tr_diff_triangle( m, idx, sign );
    for( s[m] = 0.0, n = 0; n < 4; n++ ) s[m] += sign[n] * tr[idx[n]];
  };
  for( m = 0; m < HYPERSIMPLEX_PER_SITE * param->ipw[param->D]; m++ ){
    Tr_diff_simplex( m, idx, sign );
    for( hs[m] = 0.0, n = 0; n < 5; n++ ) hs[m] += sign[n] * s[idx[n]];
    if( fabs(hs[m]) > FLT_EPSILON )
      panic("Diff triangles --> simplex --> hypersimplex FAILED");
  };
  /* ================================================================== */
  free(st); free(lnk); free(tr); free(s); free(hs);
};
/* ************************************************************************************ */
/* *******   Vertices entering various cells. Order agrees with differential  ********* */
/* ************************************************************************************ */
static void cell_constituent_poins( uint cell_idx, uint *idx, uchar dim ){
  uchar i, parity = 0;
  uint m = 0, mask = 0;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( dim == 2 ){
    if( cell_idx >= LINK_PER_SITE * param->ipw[param->D]) panic("Illegal index");
    m = cell_idx/( (uint) LINK_PER_SITE);
    get_parity( m, &parity );
    mask = l_enum_table[parity][cell_idx - LINK_PER_SITE * m];
  }else if( dim == 3 ){
    if( cell_idx >= TRIANGLE_PER_SITE * param->ipw[param->D]) panic("Illegal index");
    m = cell_idx/( (uint) TRIANGLE_PER_SITE);
    get_parity( m, &parity );
    mask = tr_enum_table[parity][cell_idx - TRIANGLE_PER_SITE * m];
  }else if( dim == 4 ){
    if( cell_idx >= SIMPLEX_PER_SITE * param->ipw[param->D]) panic("Illegal index");
    m = cell_idx/( (uint) SIMPLEX_PER_SITE );
    get_parity( m, &parity );
    mask = s_enum_table[parity][cell_idx - SIMPLEX_PER_SITE * m];
  }else if( dim == 5 ){
    if( cell_idx >= HYPERSIMPLEX_PER_SITE * param->ipw[param->D]) panic("Illegal index");
    m = cell_idx/( (uint) HYPERSIMPLEX_PER_SITE );
    get_parity( m, &parity );
    mask = reflection( hs_v_table[cell_idx - HYPERSIMPLEX_PER_SITE * m], parity, 5);
  }else
    panic("Illegal dimensionality");
  for( i = 0; i < dim; i++ )
    idx[i] = shift_index( m, (mask & (0x0F << (4*i))) >> (4*i), 1 );
};
/* ------------------------------------------------------------------------- */
void Tr_points_in_link( uint l_idx, uint *idx ){
  cell_constituent_poins( l_idx, idx, 2 );
};
/* ------------------------------------------------------------------------- */
void Tr_points_in_triangle( uint tr_idx, uint *idx ){
  cell_constituent_poins( tr_idx, idx, 3 );
};
/* ------------------------------------------------------------------------- */
void Tr_points_in_simplex( uint s_idx, uint *idx ){
  cell_constituent_poins( s_idx, idx, 4 );
};
/* ------------------------------------------------------------------------- */
void Tr_points_in_hypersimplex( uint hs_idx, uint *idx ){
  cell_constituent_poins( hs_idx, idx, 5 );
};
/* ************************************************************************************ */
/* *****  Construction of road map from sites to various cells  *********************** */
/* ************************************************************************************ */
static uint site_dependence_value( uint m_s, uchar *x, uchar i ){
  uint ret = i & INDEX_MASK;
  uchar x_s[4];
#ifdef DEBUG
  uint dis;
  uchar k, y0[4], y1[4];
#endif
  site_coordinates( x_s, m_s );
#ifdef DEBUG
  for( k = 0; k < 4; k++ ){
    bzero( y0, 4 * sizeof(uchar));
    bzero( y1, 4 * sizeof(uchar));
    y0[k] = x_s[k];
    y1[k] = x[k];
    dis = get_distance2( site_index(y0), site_index(y1) );
    if( dis && dis != 1 ) panic("Reference points are too apart");
  };
#endif
  for( i = 0; i < 4; i++ )
    if( x_s[i] != x[i] ){
      ret |= 1 << (8 + i);
      if( ( x_s[i] == 0 && x[i] != 1 ) || x_s[i] > x[i] ) ret |= 1 << (12 + i);
    };
  return ret;
};
/* ------------------------------------------------------------------------- */
static void make_site_dependence(){
  uchar j, k, parity, x[4], y[4], shift[4];
  uint m, n,idx[5], count[4];
  if( site_dependence_done ) return;
  for( m = 0; m < 4; m++ )
    for( parity = 0; parity < MAX_PARITY; parity++ ){
      site_dependence[m][parity] = (uint *) calloc( site_dependence_dim[m][parity], sizeof(uint) );
      if( !site_dependence[m][parity] ) panic("Memory");
    };
  for( parity = 0; parity < MAX_PARITY; parity++ ){
    for( m = 0; m < 4; m++ ) count[m] = 0;
    m = shift_index( 0, parity, 1 );
    site_coordinates( x, shift_index( m, 0x0F, -1) );
    for( shift[0] = 0; shift[0] < 3; shift[0]++ ){
      y[0] = (x[0] + shift[0]) % param->size[0];
      for( shift[1] = 0; shift[1] < 3; shift[1]++ ){
	y[1] = (x[1] + shift[1]) % param->size[1];
	for( shift[2] = 0; shift[2] < 3; shift[2]++ ){
	  y[2] = (x[2] + shift[2]) % param->size[2];
	  for( shift[3] = 0; shift[3] < 3; shift[3]++ ){
	    y[3] = (x[3] + shift[3]) % param->size[3];
	    n = site_index( y );
	    for( j = 0; j < LINK_PER_SITE; j++ ){
	      Tr_points_in_link( L_ENUM(j, n) , idx );
	      for( k = 0; k < 2; k++ )
		if( idx[k] == m ){
		  if( count[0] >= site_dependence_dim[0][parity] ) panic("Overflow in link counter");
		  *(site_dependence[0][parity] + count[0]) =  site_dependence_value( m, y, j );
		  count[0] += 1;
		};
	    };
	    for( j = 0; j < TRIANGLE_PER_SITE; j++ ){
	      Tr_points_in_triangle( TR_ENUM(j, n) , idx );
	      for( k = 0; k < 3; k++ )
		if( idx[k] == m ){
		  if( count[1] >= site_dependence_dim[1][parity] ) panic("Overflow in triangle counter");
		  *(site_dependence[1][parity] + count[1]) =  site_dependence_value( m, y, j );
		  count[1] += 1;
		};
	    };
	    for( j = 0; j < SIMPLEX_PER_SITE; j++ ){
	      Tr_points_in_simplex( S_ENUM(j, n) , idx );
	      for( k = 0; k < 4; k++ )
		if( idx[k] == m ){
		  if( count[2] >= site_dependence_dim[2][parity] ) panic("Overflow in simplex counter");
		  *(site_dependence[2][parity] + count[2]) =  site_dependence_value( m, y, j );
		  count[2] += 1;
		};
	    };
	    for( j = 0; j < HYPERSIMPLEX_PER_SITE; j++ ){
	      Tr_points_in_hypersimplex( HS_ENUM(j, n) , idx );
	      for( k = 0; k < 5; k++ )
		if( idx[k] == m ){
		  if( count[3] >= site_dependence_dim[3][parity] ) panic("Overflow in hypersimplex counter");
		  *(site_dependence[3][parity] + count[3]) =  site_dependence_value( m, y, j );
		  count[3] += 1;
		};
	    };
	  }; /* shift[3] */
	}; /* shift[2] */
      }; /* shift[1] */
    }; /* shift[0] */
#ifdef DEBUG
    for( m = 0; m < 4; m++ )
      if( count[m] != site_dependence_dim[m][parity] )
	panic4("Underflow in counter number %d : value %d must be %d",
	       m, count[m], site_dependence_dim[m][parity] );
#endif
  }; /* parity */
  site_dependence_done = 1;
};
/* ------------------------------------------------------------------------- */
static uint cells_dependent_on_site( uint m_s, uint *idx, uchar dim ){
  uint m, n, mask, ret = 0;
  uchar i, parity;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( m_s >= param->ipw[param->D] ) panic("Illegal parameter");
  if( !site_dependence_done ) make_site_dependence();
  get_parity( m_s, &parity );
  ret = site_dependence_dim[dim][parity];
  for( n = 0; n < ret; n++ ){
    mask = *(site_dependence[dim][parity] + n);
    for( m = m_s, i = 0; i < 4; i++ )
      if( mask & (1 << (8 + i)) )
	m = (mask & (1 << (12 + i))) ? index_down(m, i) : index_up(m, i) ;
    if( dim == 0 ){
      idx[n] = L_ENUM( mask & INDEX_MASK, m );
    }else if( dim == 1 ){
      idx[n] = TR_ENUM( mask & INDEX_MASK, m );
    }else if( dim == 2 ){
      idx[n] = S_ENUM( mask & INDEX_MASK, m );
    }else if( dim == 3 ){
      idx[n] = HS_ENUM( mask & INDEX_MASK, m );
    }else panic("Illegal parameter");
  };
  return ret;
};
/* ------------------------------------------------------------------------- */
uint Tr_link_dependent_on_site( uint m_s, uint *idx ){
  return cells_dependent_on_site( m_s, idx, 0 );
};
/* ------------------------------------------------------------------------- */
uint Tr_triangle_dependent_on_site( uint m_s, uint *idx ){
  return cells_dependent_on_site( m_s, idx, 1 );
};
/* ------------------------------------------------------------------------- */
uint Tr_simplex_dependent_on_site( uint m_s, uint *idx ){
  return cells_dependent_on_site( m_s, idx, 2 );
};
/* ------------------------------------------------------------------------- */
uint Tr_hypersimplex_dependent_on_site( uint m_s, uint *idx ){
  return cells_dependent_on_site( m_s, idx, 3 );
};
/* ************************************************************************************ */
/* *******  Restoration of normal hypercubical geometry  ****************************** */
/* ************************************************************************************ */
void Tr_hypercube_from_hypersimplex( uint m, uint *idx, int *sign ){
  char sgn;
  uchar i;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( m >= param->ipw[param->D] ) panic("Illegal parameter");
  sgn = get_parity( m, &i );
  for( i = 0; i < HYPERSIMPLEX_PER_SITE; i++ ){
    idx[i] = HS_ENUM( i, m );
    sign[i] = (i % 2) ? -sgn : sgn;
  };
};
/* ------------------------------------------------------------------------- */
void Tr_cube_from_simplex( uint m_c, uint *idx, int *sign ){
  uint m, mask, count;
  char sgn;
  uchar parity, i, i_hs, i_s, where, which , cb;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( m_c >= param->D * param->ipw[param->D] ) panic("Illegal parameter");
  m = m_c / ((uint) param->D );
  i = m_c - param->D * m;
  m = index_up( m, i );
  sgn = (-1) * get_parity( m, &parity );
  for( count = 0, i_hs = 0; i_hs < HYPERSIMPLEX_PER_SITE; i_hs++ )
    for( i_s = 0; i_s < 5; i_s++ ){
      mask = reflection( erase_bits( hs_v_table[i_hs], i_s, 5 ), parity, 4 );
      cb = common_bits( mask, &where, &which, 4 );

      if( !cb || which || where != i ) continue;

      idx[count] = s_enum( m, parity, mask);
      sign[count] = sgn;
      sign[count] *= ( i_s  %2 ) ? -1 : 1;
      sign[count] *= ( i_hs %2 ) ? -1 : 1;
      count++;
    }; /* - i_hs, i_s - */
};
/* ************************************************************************************ */
/* *****  Construction of road map from links to triangles  *************************** */
/* ************************************************************************************ */
static uint link2triangle_value( uint m_l, uchar *x, uchar i, int sign ){
  uint ret = i & INDEX_MASK;
  uchar x_l[4];
#ifdef DEBUG
  uint dis;
  uchar k, y0[4], y1[4];
#endif
  site_coordinates( x_l, m_l );
#ifdef DEBUG
  for( k = 0; k < 4; k++ ){
    bzero( y0, 4 * sizeof(uchar));
    bzero( y1, 4 * sizeof(uchar));
    y0[k] = x_l[k];
    y1[k] = x[k];
    dis = get_distance2( site_index(y0), site_index(y1) );
    if( dis && dis != 1 ) panic("Reference points are too apart");
  };
#endif
  for( i = 0; i < 4; i++ )
    if( x_l[i] != x[i] ){
      ret |= 1 << (8 + i);
      if( ( x_l[i] == 0 && x[i] != 1 ) || x_l[i] > x[i] ) ret |= 1 << (12 + i);
    };
  if( sign < 0 ) ret |= 0x0F << 16;
  return ret;
};
/* ------------------------------------------------------------------------- */
static void link2triangle_make(){
  uchar i, j, k, parity, x[4], y[4], shift[4];
  uint m, n, l, idx[3], count;
  int sign[3];
  if( link2triangle_done ) return;
  for( parity = 0; parity < MAX_PARITY; parity++ )
    for( i = 0; i < LINK_PER_SITE; i++ ){
      link2triangle[parity][i] = (uint *) calloc( link2triangle_dim[parity][i], sizeof(uint) );
      if( !link2triangle[parity][i] ) panic("Memory");
    };
  for( parity = 0; parity < MAX_PARITY; parity++ ){
    m = shift_index( 0, parity, 1 );
    site_coordinates( x, shift_index( m, 0x0F, -1) );
    for( i = 0; i < LINK_PER_SITE; i++ ){
      count = 0;
      l = L_ENUM( i, m );
      for( shift[0] = 0; shift[0] < 3; shift[0]++ ){
	y[0] = (x[0] + shift[0]) % param->size[0];
	for( shift[1] = 0; shift[1] < 3; shift[1]++ ){
	  y[1] = (x[1] + shift[1]) % param->size[1];
	  for( shift[2] = 0; shift[2] < 3; shift[2]++ ){
	    y[2] = (x[2] + shift[2]) % param->size[2];
	    for( shift[3] = 0; shift[3] < 3; shift[3]++ ){
	      y[3] = (x[3] + shift[3]) % param->size[3];
	      n = site_index( y );
	      for( j = 0; j < TRIANGLE_PER_SITE; j++ ){
		Tr_diff_link( TR_ENUM(j, n), idx, sign );
		for( k = 0; k < 3; k++ )
		if( idx[k] == l ){
		  if( count >= link2triangle_dim[parity][i] ) panic("Overflow in link2triangle counter");
		  *(link2triangle[parity][i] + count) = link2triangle_value( m, y, j, sign[k] );
		  count += 1;
		};
	      };
	    }; /* shift[3] */
	  }; /* shift[2] */
	}; /* shift[1] */
      }; /* shift[0] */
#ifdef DEBUG
      if( count != link2triangle_dim[parity][i] )
	panic3("Underflow in counter : value %d must be %d", count, link2triangle_dim[parity][i]);
#endif
    }; /* i */
  }; /* parity */
  link2triangle_done = 1;
};
/* ------------------------------------------------------------------------- */
uchar link2triangle_get( uint m_l, uint *idx, int *sign ){
  uint m, mask;
  uchar i, n, i_l, parity, ret = 0;
  if( !param ) panic("Geometry undefined");
  if( param->D != 4 ) panic("For D=4 only");
  for( i = 0; i < 4; i++ )
    if( param->size[i] % 2 ) panic("Lattice MUST be even");
  if( m_l >= LINK_PER_SITE * param->ipw[param->D] ) panic("Illegal parameter");
  if( !link2triangle_done ) link2triangle_make();
  m = m_l / ((uint) LINK_PER_SITE );
  i_l = m_l - LINK_PER_SITE * m;
  m_l = m;
  get_parity( m_l, &parity );
  ret = link2triangle_dim[parity][i_l];

  for( n = 0; n < ret; n++ ){
    mask = *(link2triangle[parity][i_l] + n);
    for( m = m_l, i = 0; i < 4; i++ )
      if( mask & (1 << (8 + i)) )
	m = (mask & (1 << (12 + i))) ? index_down(m, i) : index_up(m, i) ;
    idx[n] = TR_ENUM( mask & INDEX_MASK, m );
    sign[n] = ( mask & (0x0F << 16) ) ? -1 : 1 ;
  };
  return ret;
};
/* *********************************************************************************** */
/* ******  Self test   *************************************************************** */
/* *********************************************************************************** */
void Tr_check_all(){
  printf("Checking differential ... "); fflush(stdout);
  Tr_check_diff();
  printf("OK\n"); fflush(stdout);
};
