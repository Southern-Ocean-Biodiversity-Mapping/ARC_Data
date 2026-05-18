## load annotations, calculate total numbers and create overview file
library(SQAPI)
api <- SQAPI$new()

#### download a tally of all annotations (DOESN'T WORK YET)
## on squidle, this is the API string that shows up
## {"filters":[{"name":"point","op":"has","val":{"name":"media","op":"has","val":{"name":"deployment","op":"has","val":{"name":"platform_id","op":"in","val":[22]}}}}]}
## we change it to this in R:
filter <- query_filter(name="point",op="has",val=query_filter(name="media",op="has",val=query_filter(name="deployment",op="has",val=query_filter(name="platform_id",op="in",val=22))))
## api endpoint: /api/annotation/tally/label
params <- query_params(results_per_page = 200, limit = 200)
req <- request("GET",
               api,
               "api/annotation/tally/label",
               query_filters = filter,
               query_parameters = params,
               template = "data.csv")
dat <- parse_api(req)

#### download all annotations for a single media collection, or a single annotation set
## PS06 media collection id: 13581
## PS06 point-scoring annotations: 16742
## a single annotation set:
filter <- query_filter(name="annotation_set_id",op="eq",val=16742)
## to check if it's working
# req <- request("GET",
#                api,
#                "api/annotation",
#                query_filters = filter,
#                query_parameters = params,
#                template = "data.csv")
# dat <- parse_api(req)
req <- export(api,
               "api/annotation/export",
               query_filters = filter,
               template = "data.csv")
dat <- parse_api(req)

## all ids
## PS06-PS118, TANs, NBPs08/10/14/15, LMGs09/13, JRs, AA
## for point-scoring:
ann1.1 <- c(16742,16743, 16817, 16818, 16882, 16880, 16883)
ann1.2 <- c(16719,16738, 16741,
          16814, 16816, 18096, 16812)
ann1.3 <- c(16815, 16813,
          16820, 16819, 16745, 17547,
          17493)
ann1 <- c(ann1.1, ann1.2, ann1.3)
## for mobiles
ann2.1 <- c(17641, 17640, 17639, 17637, 17602, 17597, 17628)
ann2.2 <- c(17638, 17642, 17645, 
            17930, 17932, 18098, 17935,
            17931, 17933,
            17937, 17938, 17939, 17936,
            17589)
ann2 <- c(ann2.1, ann2.2)
## for VMEs
ann3 <- c(18077, 18078, 18086, 18087, 18601, 18600, 18074,
          18073, 18075, 18076, 
          18083, 18085, 18602, 18081,
          18084, 18082, 
          18089, 18088, 18079, 18080,
          18072)

## point annotations:
filter <- query_filter(name="annotation_set_id",op="in",val=ann1.1)
req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
dat.pt1 <- parse_api(req)
write.csv(dat.pt1, "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/SquidleASAIDAnnotations_202602_Points1.csv", row.names = FALSE)
filter <- query_filter(name="annotation_set_id",op="in",val=ann1.2)
req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
dat.pt2 <- parse_api(req)
write.csv(dat.pt2, "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/SquidleASAIDAnnotations_202602_Points2.csv", row.names = FALSE)
filter <- query_filter(name="annotation_set_id",op="in",val=ann1.3)
req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
dat.pt3 <- parse_api(req)
write.csv(dat.pt3, "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/SquidleASAIDAnnotations_202602_Points3.csv", row.names = FALSE)

## mobile annotations:
filter <- query_filter(name="annotation_set_id",op="in",val=ann2.1)
req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
dat.mob1 <- parse_api(req)
write.csv(dat.mob1, "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/SquidleASAIDAnnotations_202602_Mobiles1.csv", row.names = FALSE)
filter <- query_filter(name="annotation_set_id",op="in",val=ann2.2)
req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
dat.mob2 <- parse_api(req)
write.csv(dat.mob2, "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/SquidleASAIDAnnotations_202602_Mobiles2.csv", row.names = FALSE)

## VME annotations:
filter <- query_filter(name="annotation_set_id",op="in",val=ann3)
req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
dat.VME <- parse_api(req)
write.csv(dat.VME, "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/SquidleASAIDAnnotations_202602_VME.csv", row.names = FALSE)




########
bio.dir <- "C:/Users/jjansen/OneDrive - University of Tasmania/science/data_biological/"
dat.pt1 <- read.csv(paste0(bio.dir,"SquidleASAIDAnnotations_202602_Points1.csv"))
dat.pt2 <- read.csv(paste0(bio.dir,"SquidleASAIDAnnotations_202602_Points2.csv"))
dat.pt3 <- read.csv(paste0(bio.dir,"SquidleASAIDAnnotations_202602_Points3.csv"))
dat.pt <- rbind(dat.pt1, dat.pt2, dat.pt3)
dat.table <- table(dat.pt$label.lineage_names)[rev(order(table(dat.pt$label.lineage_names)))]
write.csv(dat.table, paste0(bio.dir,"SquidleASAIDAnnotations_202602_table_Points.csv"), row.names = FALSE)

dat.mob1 <- read.csv(paste0(bio.dir,"SquidleASAIDAnnotations_202602_Mobiles1.csv"))
dat.mob2 <- read.csv(paste0(bio.dir,"SquidleASAIDAnnotations_202602_Mobiles2.csv"))
dat.mob <- rbind(dat.mob1, dat.mob2)
dat.table.mob <- table(dat.mob$label.lineage_names)[rev(order(table(dat.mob$label.lineage_names)))]
write.csv(dat.table.mob, paste0(bio.dir,"SquidleASAIDAnnotations_202602_table_Mobiles.csv"), row.names = FALSE)

dat.vme <- read.csv(paste0(bio.dir,"SquidleASAIDAnnotations_202602_VME.csv"))
dat.table.vme <- table(dat.vme$label.lineage_names)[rev(order(table(dat.vme$label.lineage_names)))]
write.csv(dat.table.vme, paste0(bio.dir,"SquidleASAIDAnnotations_202602_table_VME.csv"), row.names = FALSE)



########################################
## download campaign specific csv files once for data paper record
# for(i in 1:length(ann1)){
#   filter <- query_filter(name="annotation_set_id",op="eq",val=ann1[i])
#   req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
#   dat.pt <- parse_api(req)
#   write.csv(dat.pt, paste0(bio.dir,"SquidleASAIDAnnotations_202602_Points_",dat.pt$point.media.deployment.campaign.name[1],".csv"), row.names = FALSE)
#   message(paste0(i," done"))
# }
# for(i in 1:length(ann2)){
#   filter <- query_filter(name="annotation_set_id",op="eq",val=ann2[i])
#   req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
#   dat.mob <- parse_api(req)
#   write.csv(dat.mob, paste0(bio.dir,"SquidleASAIDAnnotations_202602_Mobiles_",dat.mob$point.media.deployment.campaign.name[1],".csv"), row.names = FALSE)
# }
# for(i in 1:length(ann3)){
#   filter <- query_filter(name="annotation_set_id",op="eq",val=ann3[i])
#   req <- export(api, "api/annotation/export", query_filters = filter, template = "data.csv")
#   dat.vme <- parse_api(req)
#   write.csv(dat.vme, paste0(bio.dir,"SquidleASAIDAnnotations_202602_VME_",dat.vme$point.media.deployment.campaign.name[1],".csv"), row.names = FALSE)
# }

















### DOESN't WORK YET FOR CAMPAIGN OR MEDIA COLLECTION
filter <- query_filter(name="media_collection.id",op="eq",val=13581)
req <- export(api,
              "api/annotation/export",
              query_filters = filter,
              template = "data.csv")
dat <- parse_api(req)












# val = id for label to change from
filter2 <- query_filter(name = "label_id", op = "eq", val = 529)

req <- request("GET",
               api,
               "api/annotation",
               query_filters = list(filters, filter2),
               query_parameters = params,
               template = "data.csv")



{
  "filters": [
    {
      "name": "point",
      "op": "has",
      "val": {
        "name": "media",
        "op": "has",
        "val": {
          "name": "deployment",
          "op": "has",
          "val": {
            "name": "platform_id",
            "op": "in",
            "val": [
              22
            ]
          }
        }
      }
    }
  ]
}



