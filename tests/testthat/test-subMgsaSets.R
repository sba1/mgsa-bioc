test_that("subMgsaSets() gives correct result", {
  test.sets <- list(
      I = c( "a", "b", "c", "d", "e" ),
      II = c( "a", "b" ),
      III = c( "a", "b", "d" ),
      IV = c( "a", "c", "e" ),
      V = c( "e" ),
      VI = c( "b", "e" ) )
  test.item_anno <- data.frame( x = 1:4, row.names = c( "d", "b", "c", "a" ) )
  test.set_anno <- data.frame( y = 1:4, row.names = c( "III", "V", "I", "VI" ) )
  
  mgsa.simple <- new( "MgsaSets", sets=test.sets,
                      itemAnnotations = test.item_anno,
                      setAnnotations = test.set_anno )
  mgsa.sub <- subMgsaSets( mgsa.simple, c( "a", "c", "d" ) )
  
  test.sub_sets <- list(
      I = c( "a", "c", "d" ),
      II = c( "a" ),
      III = c( "a", "d" ),
      IV = c( "a", "c" )
  )
  expect_equal( length(mgsa.sub), length(test.sub_sets) )
  expect_equal( mgsa.sub@sets, lapply(test.sub_sets,function(items) as.integer(mgsa.sub@itemName2ItemIndex[items])) )
  expect_equal( mgsa.sub@itemName2ItemIndex, c( a = 1, c = 2, d = 3 ) )
  expect_equal( as.integer(itemIndices(mgsa.sub, c("a","b","c","d","e" )) ), c( 1, NA, 2, 3, NA ) )
  expect_equal( nrow(itemAnnotations(mgsa.sub)), 3 )
  expect_equal( itemAnnotations(mgsa.sub, c("a","b","c","d","e" ) ),
                data.frame( x = c( 4, NA, 3, 1, NA ),
                            row.names = c( 'a', 'NA', 'c', 'd', 'NA.1' ) ) )
  expect_equal( nrow(setAnnotations(mgsa.sub)), 2 )
  expect_equal( setAnnotations(mgsa.sub, c("I","II","III","IV","V", "VI" ) ),
                data.frame( y = c( 3, NA, 1, NA, NA, NA ),
                            row.names = c( 'I', 'NA', 'III', 'NA.1', 'NA.2', 'NA.3' ) ) )
} )