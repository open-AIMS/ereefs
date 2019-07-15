#' Download river discharge (flow) data from Qld. Govt. servers at https://water-monitoring.information.qld.gov.au
#'
#' @param station Guaging station ID (obtainable from https://water-monitoring.information.qld.gov.au). Default is "124001B".
#'        Stations used in eReefs models are:
#'                oconnell: "124001B")
#'                tully: "113006A"
#'                pioneer: "125016A"
#'                burnett: "136007A"
#'                fitzroy: "130005A"
#'                herbert: "116001F"
#'                normanby: "105107A"
#'                burdekin: "120006B"
#'                mulgrave: "111007A"
#'                russell: "111101D"
#'                barron: "110001D"
#'                daintree: "108002A"
#'                sjohnstone: "112101B"
#'                njohnstone: "112004A"
#' @param start_time Date from which to extract data, as text "YYYYMMDD". Default is "20101001"
#' @param end_time Last date required, as text "YYYYMMDD". Default is "20180930"
#' @param interval Time interval over which to average discharge, e.g. "year", "week", "day". Default is "year". 
#' @return  mean discharge in ML/d
#' @export
download_discharge <- function(station="124001B", start_time="20101001", end_time="20180930", interval="year") {
   params <- paste0("{params{", 
                    "&site_list=", station, 
                    "&start_time=", start_time, "000000",
                    "&varfrom=140.00",
                    "&interval=", interval,
                    "&varto=141.00",
                    "&datasource=AT",
                    "&end_time=", end_time, "000000",
                    "&data_type=mean",
                    "&multiplier=1",
                    "&{function=get_ts_traces",
                    "&version=2}}")

   wmip_data <- jsonlite::fromJSON(txt=paste0("https://water-monitoring.information.qld.gov.au/cgi/webservice.pl?", params))
   discharge <- as.numeric(wmip_data$return['traces']$trace$trace[[1]]$v)
   quality <- wmip_data$return['traces']$trace$trace[[1]]$q
   datevec <- as.Date(strptime(wmip_data$return['traces']$trace$trace[[1]]$t / 1e6, format="%Y%m%d"))
   return(data.frame(date = datevec, discharge=discharge, quality=quality))
}
