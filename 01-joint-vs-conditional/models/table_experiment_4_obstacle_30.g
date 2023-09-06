Include: 'table_experiment_4_obstacle.g'
Edit goal_table{ size:[.30,.30, .04, .01] }

obs_2(goal_table) {
joint:rigid, shape:ssBox, Q:[.1 .15 .065 1 0 0 1], size:[.1 .1 .09 .01], color:[.6 .6 .6 1] 
}

obs_3(goal_table) {
joint:rigid, shape:ssBox, Q:[-.1 .2 .065 1 0 0 1], size:[.1 .1 .09 .01], color:[.6 .6 .6 1] 
}


