#!/usr/bin/env anduril
//$OPT --threads 75
//$OPT --pipe "tee $ANDURIL_EXECUTION_DIR/_log"
//$OPT --pipe "anduril-pager --ls"
//$OPT --wrapper slurm-prefix
//$PRE export ANDURIL_SELECTNODE=plain

import anduril.builtin._
import anduril.tools._
import org.anduril.runtime._
import anduril.microarray._
import anduril.sequencing._

object fits_kmeans{

    val k = 25
    val n = 10

    val seeds = REvaluate(param1 = k.toString, 
                          param2 = n.toString, 
                          script = """
                                    runs = rep(seq.int(as.numeric(param1)), each = as.numeric(param2))
                                    set.seed(901012)
                                    seeds = sample.int(5000, length(runs))
                                    table.out = data.frame(k = runs, seed = seeds, key = seq.int(length(runs) )) 
                                   """)   
    
    val script = StringInput(content="""
                 source("poi_decom_gain/poi_decom.R", chdir = T)
                 weights <- readRDS("weights.RDS")
                 gains <- readRDS("gains.RDS")
                 counts <- readRDS("counts.RDS")
                 
                 set.seed(as.numeric(param2))
                 fit <- poi_decom_control(Y = counts, D = weights, g = gains, k = as.numeric(param1), max.iter = 200)

                 table.out = data.frame()
                 setwd(document.dir)
                 saveRDS(fit, file = paste0("fit_", param1, ".RDS"))
                 
        """)    

    val fits  = NamedMap[REvaluate]("fits") 

    for (row <- iterCSV(seeds.table)) {	  
        val k = row("k")             
        val seed = row("seed")
        val key = row("key")
        fits(key) = REvaluate(param1 = k,                          
                              param2 = seed,             
                              script = script)  

   }

}
