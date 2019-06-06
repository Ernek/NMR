%chk=test_run                    ! Checkpoint file needed to save info needed to restart job    
#P hf/cc-pVDZ                    ! Job directives . Hartree Fock single point calc with cc-pVDZ basis sets
                                 ! Empty line ... This empty line is needed
job name                         ! Job name descriptor
                                 ! Empty line ... This empty line is needed
-1 1                             ! Charge Spin
Al 4.158934 12.319204 8.139289   ! Atomic-symbol x-coord y-coord z-coord
O 5.204656 12.33906 6.772874
H 5.524529 11.426276 6.453445
O 4.893735 13.419199 9.302062
H 5.82317 13.510635 9.461939
O 2.621193 13.123312 7.624101
H 2.124667 13.404712 8.410292
O 3.723529 10.689556 8.948275
H 3.223616 10.091097 8.387248
                                 ! Empty line ... This empty line is needed  
