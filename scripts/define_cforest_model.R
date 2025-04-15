
library(parsnip)
library(rlang)
library(partykit)

cforest_spec <- function(mtry = NULL, minsplit = NULL, ntree = NULL,formula=NULL ,...) {
  parsnip::new_model_spec(
    "cforest",  # Model name
    args = list(
      mtry = rlang::enquo(mtry),
      minsplit = rlang::enquo(minsplit),
      ntree = ntree,
      formula=formula,
      ... = rlang::enquos(...)
    ),
    method = "cforest",  # Underlying method
    mode = "regression",  # Regression mode
    eng_args = list(
      control = ctree_control(
        teststat = "quad", 
        testtype = "Univ", 
        mincriterion = 0, 
        saveinfo = FALSE,
        ...
      ),
      mtry = mtry,  # Pass additional arguments for engine
      ntree = ntree
    )  # Engine-specific arguments
  )
}

set_new_model("cforest")
set_model_mode("cforest", mode = "regression")
set_model_engine("cforest", mode = "regression", eng = "partykit")
set_dependency("cforest", eng = "partykit", pkg = "partykit")

set_fit(
  model = "cforest",
  eng = "partykit",
  mode = "regression",
  value = list(
    interface = "formula",
    protect = c("formula", "data"),
    func = c(pkg = "partykit", fun = "cforest"),
    defaults = list()
  )
)

set_encoding(
  model = "cforest",
  eng = "partykit",
  mode = "regression",
  options = list(
    predictor_indicators = "none",
    compute_intercept = FALSE,
    remove_intercept = FALSE,
    allow_sparse_x = FALSE
  )
)

