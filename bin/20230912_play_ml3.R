
library(mlr3)

task = tsk("penguins")
split = partition(task)
learner = lrn("classif.rpart")

learner$train(task, row_ids = split$train)
learner$model

prediction = learner$predict(task, row_ids = split$test)
prediction

prediction$score(msr("classif.acc"))

#

library(mlr3verse)

tasks = tsks(c("breast_cancer", "sonar"))

glrn_rf_tuned = as_learner(
  ppl("robustify") %>%
    auto_tuner(
      tnr("grid_search", resolution = 5),
      lrn("classif.ranger", num.trees = to_tune(200, 500)),
      rsmp("holdout")
    )
)
glrn_rf_tuned$id = "RF"
