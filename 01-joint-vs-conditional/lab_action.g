

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

