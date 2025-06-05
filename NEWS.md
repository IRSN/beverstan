
**beverstan** Package News
======================

# New in version 0.0.3

## Enhancements

- The `"TSGEVBayes"` class has now a `profLik` method that can be used
  to profile one or several function(s) of the parameter such as a
  quantile or a return period. Due to the re-definition of the S3
  generic function `profLik`, a conflict with the **`{NSGEV}`** package
  exists. So for now the methods must be called by their names as
  `beverstan:::profLik.TVGEVBayes`. This will be fixed in a future
  version.
  
- Some methods such as `autoplot` have been implemented for the class
  `"profLik.TVGEVBayes"` of the objects returned by the `profLik`
  method.
  
     
