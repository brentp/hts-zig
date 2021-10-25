[![Build](https://github.com/brentp/hts-zig/actions/workflows/build.yml/badge.svg)](https://github.com/brentp/hts-zig/actions/workflows/build.yml)

# WIP [zig](https://ziglang.org) wrapper for [htslib](htslib.org) parsing of VCFs for genetic variants

I wrote this learning zig so it probably has many non-ziggy usages.
Much of htslib VCF stuff is *not* supported, but some of the difficult parts are
implemented the rest is fairly mechanical.


# Usage

```zig
const VCF = @import("src/vcf.zig").VCF
const allocator = std.testing.allocator;

var vcf = VCF.open("tests/test.snpeff.bcf").?; // can return null
var variant = ivcf.next().?;
try stdout.print("\nvariant:{any}\n", .{variant}); // Variant(chr1:30859-30860 (G/C))

// get the AD field
var fld = "AD";

// needs to allocate, this interface will likely change.

// extract the FORMAT/sample AD field (allelic depth)
var ad = try variant.samples(i32, fld, allocator);
// 4 samples * 2
try std.testing.expect(ad.len == 8);
try stdout.print("\nAD:{any}\n", .{ad}); // { 7, 0, 2, 0, 6, 0, 4, 0 }

// free the memory.
allocator.free(ad);
```

# testing and dev

This will require `zig` and `htslib` installed.
Then the startig point is:
```
zig build test
```

# TODO (maybe)

+ Add `vcf.query()` (currently only iteration from start of file is supported,
  not querying by genomic location.
+ writing. currently everything is read-only.
+ other parts of htslib like bam/cram support.
+ fix ergonomics and think about error and null return types.

# Why?

I find that it's quite useful to learn a new language by writing something that
I'll use every day. I think zig looks interesting and sane
