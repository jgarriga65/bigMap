
# ------------------------------------------------------------------------------
# source xtrs
# ------------------------------------------------------------------------------

bdm.source <- function(file)
{
	ver <- getNamespaceVersion('bigMap')
	ver <- paste(substr(ver, 1, 4), as.numeric(substr(ver, 5, 5)) -1, sep = '')
	source(paste('~/bigMap/bigMap_', ver, '/R/bdm_auxf.R', sep = ''))
	source(paste('~/bigMap/xtrs/', file, '.R', sep = ''))
}
