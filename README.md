[![Build](https://github.com/brentp/hts-zig/actions/workflows/build.yml/badge.svg)](https://github.com/brentp/hts-zig/actions/workflows/build.yml)

# WIP [zig](https://ziglang.org) wrapper for [htslib](htslib.org) parsing of VCFs for genetic variants

I wrote this learning zig so it probably has many non-ziggy usages.
Most of the htslib VCF/Variant stuff is supported.

# Usage

> :warning: `hts-zig` tries to allocate as little as possible, use
> `.dup()` on a variant if when it's necessary to keep it in memory.

```zig
const vcf = @import("src/vcf.zig")
const VCF = vcf.VCF

const allocator = std.testing.allocator;

var vcf = VCF.open("tests/test.snpeff.bcf").?; // can return null
// this will close the file and cleanup htslib memory at end of scope
defer vcf.deinit(); 
var variant = vcf.next().?; # .? syntax gets an optional result.
try stdout.print("\nvariant:{any}\n", .{variant}); // Variant(chr1:30859-30860 (G/C))

// # get the AD field
// # needs to allocate, this interface will likely change.
// # extract the FORMAT/sample AD field (allelic depth)
// # to get INFO, use get(vcf.Field.info, ...);
var fld = "AD";
var ads = std.ArrayList(i32).init(allocator)
defer ads.deinit();
try variant.get(vcf.Field.format, ads, fld);
// 4 samples * 2
try stdout.print("\nAD:{any}\n", .{ads.items}); // { 7, 0, 2, 0, 6, 0, 4, 0 }

// # genotypes:
var gts_mem = std.ArrayList(i32).init(allocator); // can re-use this for each variant
defer gts_mem.deinit();
var gts = try variant.genotypes(&gts_mem);
try stdout.print("\ngts:{any}\n", .{gts});
// # gts:[0/0/, 0/0/, 0/1/, 0/0/] (note trailing / is accurate as it's how it's stored in htslib)

// # region queries
var iter = try ivcf.query(chrom, 69269, 69270);
while (iter.next()) |v| {
    try std.testing.expect(v.start() == 69269);
}
```

# testing and dev

This will require `zig` and `htslib` installed.
```
zig build test
```
Zig requires very little additional work to wrap a C library, but it is
currently not able to access struct fields when the struct contains bitfields.
For this reason, we need to write functions to access fields in many htslib
structs. Those functions are [here](https://github.com/brentp/hts-zig/blob/main/src/htslib_struct_access.c)

# TODO

- [X] Add nice genotypes access methods/structs
- [X] Add `vcf.query()` (currently only iteration from start of file is supported, not querying by genomic location.
- [X] updating header.
- [X] setting INFO fields.
- [X] writing. currently everything is read-only.
- [X] fewer allocations (pass ArrayList to functions).
- [X] support querying vcf as well as bcf
- [X] set_samples().
- [ ] fix ergonomics and think about error and null return types.

# Why?

it's quite useful to learn a new language by writing something that
I'll use every day. zig looks interesting and sane
