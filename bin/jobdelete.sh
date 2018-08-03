#!/bin/bash 

             COUNTER=1287138
                      while [  $COUNTER -lt 1287292 ]; do
                                   qdel $COUNTER
                                                let COUNTER=COUNTER+1 
                                                         done