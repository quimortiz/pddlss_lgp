Include: 'table_experiment.g'
Edit goal_table { size: [0.1, 0.1, 0.04 , 0.01 ] }

obs_1(goal_table) {
joint:rigid, shape:ssBox, Q:[.1 -0.05 .065 1 0 0 1], size:[.1 .1 .09 .01], color:[.6 .6 .6 1] 
}

obs_2(goal_table) {
joint:rigid, shape:ssBox, Q:[.1 .1 .065 1 0 0 1], size:[.1 .1 .09 .01], color:[.6 .6 .6 1] 
}

obs_3(goal_table) {
joint:rigid, shape:ssBox, Q:[-.12 -.1 .065 1 0 0 1], size:[.1 .1 .09 .01], color:[.6 .6 .6 1] 
}


