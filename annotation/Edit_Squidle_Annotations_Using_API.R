library(remotes)
install_github("sajessop/SQAPI")
library(SQAPI)
api <- SQAPI$new()


## GET request to check first
# val = id for annotation set
filters <- query_filter(name = "annotation_set_id", op = "eq", val = 17935) ## PS81 is 17602, PS14 is 17640, PS18 is 17639, PS96 is 17597, PS118 is 17628, AA2011 is 17589
# val = id for label to change from
filter2 <- query_filter(name = "label_id", op = "eq", val = 15005)
params <- query_params(results_per_page = 200, limit = 200)
req <- request("GET",
               api,
               "api/annotation",
               query_filters = list(filters, filter2),
               query_parameters = params,
               template = "data.csv")
dat <- parse_api(req)


## PATCH request
# JSON PATCH DATA:
### ID:529 is for brittle / snake star
patch <- list("label_id" = 529)

for(i in 1:10){
  message(i)
  patch_req <- request(
  "PATCH",
  api,
  "api/annotation",
  query_filters = list(filters, filter2),
  query_parameters = params,
  body = patch
)
}
