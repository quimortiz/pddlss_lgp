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




poseEq
touch
quimgrasp
impulse
stable
stableOn
dynamic
dynamicOn
liftDownUp


# objects
block1 , block2 
block3
block1_ref,  block2_ref, block3_ref
r_gripper 
l_gripper 
goal1_table
# goal2_table
goal3_table
# obstacle


## initial state (generated by the code)
START_STATE {


(is_box goal1_table)
(is_sky goal1_table) 



(is_box goal3_table)
(is_sky goal3_table) 

# (is_box goal2_table)
# (is_sky goal2_table) 

# whether the robot can move or not the big boxes.

(is_object block1) 
(is_object block2) 
(is_object block3) 



# (is_box block1) 
# (is_box block2) 


(is_sky block1) 
(is_sky block2) 
(is_sky block3) 





# (is_box block1_ref) 
# (is_box block2_ref) 
# (is_box block3_ref) 
# (is_box block4_ref) 
# (is_box block5_ref) 
# (is_box block6_ref) 

# new obstacle
# (is_box obs)
# (is_box obs_ref) 



(is_gripper r_gripper)
(is_gripper l_gripper)


(stable block1_ref block1 )
(stable block2_ref block2 )
(stable block3_ref block3 )



(on block1_ref block1 )
(on block2_ref block2 )
(on block3_ref block3 )


(poseEq block1_ref block1 )
(poseEq block2_ref block2 )
(poseEq block3_ref block3 )






}



### RULES

#####################################################################

### Reward
REWARD {
}

#####################################################################


# Include: 'lab_action.g'
Include: 'lab_action_quim.g'


#####################################################################

