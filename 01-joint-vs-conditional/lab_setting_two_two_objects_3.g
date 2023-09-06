# Include: '../../rai-robotModels/scenarios/pandasTable-calibrated-onlyone.g'
Include: '../rai-robotModels/scenarios/pandasTable-calibrated.g'

right_marker(r_gripper) { shape:sphere, size:[.2], color:[1,0,0,.5]}


goal1_table (table){ joint:rigid, shape:ssBox, Q:[.0 .2 .035 1 0 0 0], size:[.5 .5 .04 .01], color:[.1 .9 .1] }
# goal2_table (table){ joint:rigid, shape:ssBox, Q:[1 .2 .035 1 0 0 0], size:[.2 .2 .04 .01], color:[.9 .1 .1] }
goal3_table (table){ joint:rigid, shape:ssBox, Q:[-1 .2 .035 1 0 0 0], size:[.3 .3 .04 .01], color:[.1 .1 .9] }


block1_ref (table){ joint:rigid, shape:ssBox, Q:[.8 0.0 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.6 .6 .8 .2] }
block2_ref (table){ joint:rigid, shape:ssBox, Q:[.8 .2 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.8 .6 .8 .2] }
block3_ref (table){ joint:rigid, shape:ssBox, Q:[.8 .4 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.8 .6 .6 .2] }
# block4_ref (table){ joint:rigid, shape:ssBox, Q:[-.8 .1 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.6 .8 .6 .2] }
# block5_ref (table){ joint:rigid, shape:ssBox, Q:[-.8 .3 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.6 .8 .8 .2] }
# block6_ref (table){ joint:rigid, shape:ssBox, Q:[-.8 .5 .095 1 0 0 1], size:[.06 .15 .09 .01], color:[.8 .8 .8 .2] }

# Blocks
# NOTE: should this depend on 
block1(block1_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.6 .6 .8] }
block2(block2_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.8 .6 .8] }
block3(block3_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.8 .6 .6] }
# block4(block4_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.6 .8 .6] } 
# block5(block5_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.6 .8 .8] }
# block6(block6_ref){ joint:free, shape:ssBox, size:[.06 .15 .09 .01], color:[.8 .8 .8] }

fgoal(table){ Q:<t(.5 .5 .095) d(90 0 0 1) d(90 0 1 0) > shape:marker,size:[.1]}
fgoal2(table){ Q:<t(-.4 .2 .095) d( 180 1 0  0) d(180 0 0 1)>] shape:marker,size:[.1]}

f1(block1){ shape:marker,size:[.1]}
f2(block2){ shape:marker,size:[.1]}
f3(block3){ shape:marker,size:[.1]}
# f4(block4){ shape:marker,size:[.1]}
# f5(block5){ shape:marker,size:[.1]}
# f6(block6){ shape:marker,size:[.1]}

Delete camera
