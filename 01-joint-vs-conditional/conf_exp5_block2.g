# Include: '../../rai-robotModels/scenarios/pandasTable-calibrated-onlyone.g'
# Include: '../rai-robotModels/scenarios/pandasTable-calibrated.g'

Include: '../rai-robotModels/scenarios/pandasTable-calibrated-onlyone.g'

Edit r_gripper{logical:{gripper}}



other_table (table){ joint:rigid, shape:ssBox, Q:[0 .2 .035 1 0 0 0], size:[.15 .15 .04 .01], color:[.9 .1 .1] }

goal_table (table){ joint:rigid, shape:ssBox, Q:[0.4 .2 .035 1 0 0 0], size:[.15 .15 .04 .01], color:[.9 .1 .1] }

block1_ref (table){ joint:rigid, shape:ssBox, Q:[.8 0. .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.6 .6 .8 .2] }
block2_ref (table){ joint:rigid, shape:ssBox, Q:[.4 .2 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.8 .6 .8 .2] }

block3_ref (table){ joint:rigid, shape:ssBox, Q:[.8 .2 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.8 .6 .8 .2] }


block1(block1_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.6 .6 .8] }
block2(block2_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.8 .6 .8] } block3(block3_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.6 .8 .6] }

f1(block1){ shape:marker,size:[.1]}
f2(block2){ shape:marker,size:[.1]}
f3(block3){ shape:marker,size:[.1]}

f1_rand(goal_table) { shape:marker , Q:<t(0.01 0 0.1) d(-0.5 0 0 1)> , size:[.1] }
f2_rand(goal_table) { shape:marker , Q:<t(0.0 0.01 0.1) d(0.5 0 0 1)> ,  size:[.1] }
f3_rand(goal_table) { shape:marker , Q:<t(0.0 0.01 0.1) d(0.5 0 0 1)> ,  size:[.1] }

f_block1_col(block1){  shape:ssBox, size:[.1 .19 .09 .01], color:[.6 .6 .8, .5] }
f_block2_col(block2){  shape:ssBox, size:[.1 .19 .09 .01], color:[.8 .6 .8, .5] }
f_block3_col(block3){  shape:ssBox, size:[.1 .19 .09 .01], color:[.6 .8 .6, .5] }


Delete camera

