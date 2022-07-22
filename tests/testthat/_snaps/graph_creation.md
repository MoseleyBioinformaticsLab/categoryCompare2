# jaccard graph is correct

    Code
      test_jaccard@edgeData
    Output
      An object of class "attrData"
      Slot "data":
      $`g1|g2`
      $`g1|g2`$weight
      [1] 1
      
      
      $`g1|g3`
      $`g1|g3`$weight
      [1] 0.5
      
      
      $`g2|g3`
      $`g2|g3`$weight
      [1] 0.5
      
      
      $`g1|g4`
      $`g1|g4`$weight
      [1] 0.25
      
      
      $`g2|g4`
      $`g2|g4`$weight
      [1] 0.25
      
      
      $`g3|g4`
      $`g3|g4`$weight
      [1] 0.4
      
      
      $`g4|g5`
      $`g4|g5`$weight
      [1] 0.1666667
      
      
      $`g5|g6`
      $`g5|g6`$weight
      [1] 0.6666667
      
      
      
      Slot "defaults":
      $weight
      [1] 0
      attr(,"class")
      [1] "FLOATING"
      
      

# overlap graph is correct

    Code
      test_overlap@edgeData
    Output
      An object of class "attrData"
      Slot "data":
      $`g1|g2`
      $`g1|g2`$weight
      [1] 1
      
      
      $`g1|g3`
      $`g1|g3`$weight
      [1] 1
      
      
      $`g2|g3`
      $`g2|g3`$weight
      [1] 1
      
      
      $`g1|g4`
      $`g1|g4`$weight
      [1] 0.5
      
      
      $`g2|g4`
      $`g2|g4`$weight
      [1] 0.5
      
      
      $`g3|g4`
      $`g3|g4`$weight
      [1] 0.6666667
      
      
      $`g4|g5`
      $`g4|g5`$weight
      [1] 0.3333333
      
      
      $`g5|g6`
      $`g5|g6`$weight
      [1] 1
      
      
      
      Slot "defaults":
      $weight
      [1] 0
      attr(,"class")
      [1] "FLOATING"
      
      

# combined graph is correct

    Code
      test_combined@edgeData
    Output
      An object of class "attrData"
      Slot "data":
      $`g1|g2`
      $`g1|g2`$weight
      [1] 1
      
      
      $`g1|g3`
      $`g1|g3`$weight
      [1] 0.75
      
      
      $`g2|g3`
      $`g2|g3`$weight
      [1] 0.75
      
      
      $`g1|g4`
      $`g1|g4`$weight
      [1] 0.375
      
      
      $`g2|g4`
      $`g2|g4`$weight
      [1] 0.375
      
      
      $`g3|g4`
      $`g3|g4`$weight
      [1] 0.5333333
      
      
      $`g4|g5`
      $`g4|g5`$weight
      [1] 0.25
      
      
      $`g5|g6`
      $`g5|g6`$weight
      [1] 0.8333333
      
      
      
      Slot "defaults":
      $weight
      [1] 0
      attr(,"class")
      [1] "FLOATING"
      
      

