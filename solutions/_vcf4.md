the number of variant  can be computed using this simple command:

```
grep -v "^#" variants/NA12878.rmdup.realign.hc.vcf | wc -l 
grep -v "^#" variants/NA12878.hc.vcf | wc -l 
```

We have 404 variants in the realigned vcf and 405 in the raw vcf

In that case the impact of recalibration is very low because HaplotypeCaller perform a similar step internally


