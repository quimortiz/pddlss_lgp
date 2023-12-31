QUIT
WAIT
ANY
Terminate

FOL_World{
  hasWait=false
  gamma = 1.
  stepCost = 1.
  timeCost = 0.
}

## basic predicates
is_ref
is_gripper
is_object
is_box
is_sky # it does not have anything on top
is_picked

on
busy     # involved in an ongoing (durative) activity

INFEASIBLE
INFEASIBLE_pick

## KOMO symbols
above


fgoal
fgoal2


poseEq
touch
impulse
stable
stableOn
dynamic
dynamicOn
liftDownUp


# objects
table, table_ref # i have added table. is it fine?
block1 , block2 , block3 , block4, block5 , block6
block1_ref,  block2_ref, block3_ref , block4_ref, block5_ref,  block6_ref 
l_gripper, r_gripper 
goal1_table
goal2_table
goal3_table
# obstacle


## initial state (generated by the code)
START_STATE {

(is_ref fgoal)
(is_ref fgoal2)
(is_box table)

(is_box goal1_table)
(is_box goal2_table)
(is_box goal3_table)


# whether the robot can move or not the big boxes.

(is_object block1) 
(is_object block2) 
(is_object block3) 
(is_object block4) 
(is_object block5) 
(is_object block6) 

# (is_box table_ref)


(is_box block1) 
(is_box block2) 
(is_box block3) 
(is_box block4) 
(is_box block5) 
(is_box block6) 


(is_sky block1) 
(is_sky block2) 
(is_sky block3) 
(is_sky block4) 
(is_sky block5) 
(is_sky block6) 





# (is_box block1_ref) 
# (is_box block2_ref) 
# (is_box block3_ref) 
# (is_box block4_ref) 
# (is_box block5_ref) 
# (is_box block6_ref) 

# new obstacle
# (is_box obs)
# (is_box obs_ref) 



(is_gripper l_gripper) 
(is_gripper r_gripper)


(stable block1_ref block1 )
(stable block2_ref block2 )

(stable block3_ref block3 )
(stable block4_ref block4 )
(stable block5_ref block5 )
(stable block6_ref block6 )


(on table_ref table )
(on block1_ref block1 )
(on block2_ref block2 )

(on block3_ref block3 )
(on block4_ref block4 )
(on block5_ref block5 )
(on block6_ref block6 )

(poseEq block1_ref block1 )
(poseEq block2_ref block2 )

(poseEq block3_ref block3 )
(poseEq block4_ref block4 )
(poseEq block5_ref block5 )
(poseEq block6_ref block6 )

(is_sky table)
(is_sky goal1_table)
(is_sky goal2_table)
(is_sky goal3_table)




}



### RULES

#####################################################################

### Reward
REWARD {
}

#####################################################################

DecisionRule pick {
  Obj, From, Hand, 
  { (is_gripper Hand) (is_object Obj) 
  (is_sky Obj) 
  (on From Obj) (busy Hand)! (INFEASIBLE_pick Hand Obj)! }
  { (on From Obj)! (above Obj From)! (touch From Obj)! (stable From Obj)! (stableOn From Obj)!
    (on Hand Obj)
    (poseEq From Obj)!
    (busy Hand) #logic
    (touch Hand Obj) (stable Hand Obj)  #geometric
    (is_sky From)
    # (is_sky Obj)! #NOTE: this block the handover
    (is_picked Obj)
    }
}

#####################################################################

DecisionRule place {
  Obj, Hand, To,
  { (is_gripper Hand) (is_object Obj) (on Hand Obj) (is_box To) (is_sky To) (is_picked To)! }
  { (busy Hand)! (on Hand Obj)! (stable Hand Obj)! (touch Hand Obj)!
    (on To Obj) #logic
    (above Obj To) (stableOn To Obj) #geometric
    (is_sky To)!
    (is_sky Obj)
    (INFEASIBLE_pick ANY Obj)! block(INFEASIBLE_pick ANY Obj)
    (is_picked Obj)!
    }
}


DecisionRule placeeq {
  Obj, Hand, To,
  { (is_gripper Hand) (is_object Obj) (on Hand Obj) (is_ref To)  (is_picked To)! }
  { (busy Hand)! (on Hand Obj)! (stable Hand Obj)! (touch Hand Obj)!
    (on To Obj) #logic
    # (above Obj To) (stableOn To Obj) #geometric
    (poseEq Obj To) (stable To Obj) #geometric
    (is_sky To)!
    (is_sky Obj)
    (INFEASIBLE_pick ANY Obj)! block(INFEASIBLE_pick ANY Obj)
    (is_picked Obj)!
    }
}




#####################################################################

