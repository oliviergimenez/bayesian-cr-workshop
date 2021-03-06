library(rjson) #import and manipulate JSON files
library(dplyr) #data handling / pipe char
library(kableExtra)

# Code from https://github.com/omgitsjessie/slack-exports adapted to format a Q&A

#Unzip the slack export file
#Put the unzipped folder in your R working directory (getwd() if you aren't sure where that is)
#Copy the name of the unzipped file, change the next variable to match.
#exportname <- "exportunzip" 
exportname <- "json-export"
working_directory <- getwd() %>% as.character()
slackexport_folder_path <- paste0(working_directory,"/",exportname)
slackexport_folder_path

#Make a list of all channels present in the slack export
#This information is all in the "<path>/exportname/channels.json" file
channels_path <- paste0(slackexport_folder_path,"/channels.json")
channels_json <- fromJSON(file = channels_path)
channel_list <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                          c("ch_id", "name", "created", "creator", "is_archived",
                            "is_general", "members", "topic", "purpose"))


#Create channel_list df with channel information, and add a file list into channels_json
for (channel in 1:length(channels_json)) { 
  #Make a df (channel_list) with information about each channel, from the JSON file
  channel_list[channel, "ch_id"] <- channels_json[[channel]]$id
  channel_list[channel, "name"] <- channels_json[[channel]]$name
  channel_list[channel, "created"] <- channels_json[[channel]]$created
  channel_list[channel, "creator"] <- channels_json[[channel]]$creator
  channel_list[channel, "is_archived"] <- channels_json[[channel]]$is_archived
  channel_list[channel, "is_general"] <- channels_json[[channel]]$is_general
  #make a comma separated list of members
    memberlist <- ""
    for(member in 1:length(channels_json[[channel]]$members)) {
      #if it isn't the last member
      if(member < length(channels_json[[channel]]$members)) {
        memberlist <- paste0(memberlist, channels_json[[channel]]$members[[member]], ", ")
      }
      if(member == length(channels_json[[channel]]$members)) {
        memberlist <- paste0(memberlist, channels_json[[channel]]$members[[member]])
      }
      
    }
  channel_list[channel, "members"] <- memberlist
  channel_list[channel, "topic"] <- channels_json[[channel]]$topic$value
  channel_list[channel, "purpose"] <- channels_json[[channel]]$purpose$value

  #For each channel make a list of all the individual JSON files (one file per day of activity)
  #Add that list to the channels_json object as a list in each channel: channels_json[[channel]]$dayslist
  channel_folder_path <- ""
  channels_json[[channel]]$dayslist <- ""
  channel_folder_path <- paste0(slackexport_folder_path,"/",channel_list[channel,"name"])
  channels_json[[channel]]$dayslist <- list.files(channel_folder_path, 
                                                  pattern=NULL, all.files=FALSE, full.names=FALSE)
  
}


#function to convert a single JSON file into a dataframe with specific fields extracted
slack_json_to_dataframe <- function(slack_json) {
  #blank table with correct colnames:
  messages_df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), 
                            c("msg_id", "ts", "user", "type", "text", "reply_count",
                              "reply_users_count", "ts_latest_reply", "ts_thread", 
                              "parent_user_id"))
  #for each slack message (list item in JSON file), extract relevant fields 
  for (message in 1:length(slack_json)) { 
    #messages with a file attached have no msg_id, just a file ID. Grab file ID if they have it, otherwise msg ID.
    #TODO - something is wrong with messages_df$msg_id - it is recording as NA for all obs.
    if (is.null(slack_json[[message]]$files$id) == FALSE) {
      messages_df[message, "msg_id"] <- slack_json[[message]]$files$id
    }
    if (is.null(slack_json[[message]]$msg_id) == FALSE) { 
      messages_df[message, "msg_id"] <- slack_json[[message]]$client_msg_id
    } 
    messages_df[message, "ts"] <- slack_json[[message]]$ts
    messages_df[message, "user"] <- slack_json[[message]]$user
    messages_df[message, "type"] <- slack_json[[message]]$type
    messages_df[message, "text"] <- slack_json[[message]]$text
    #Some values only occur for parents or children of threads.
    #this will trigger for all parent messages
    if (is.null(slack_json[[message]]$reply_count) == FALSE) { 
      messages_df[message, "reply_count"] <- slack_json[[message]]$reply_count
      messages_df[message, "reply_users_count"] <- slack_json[[message]]$reply_users_count
      messages_df[message, "ts_latest_reply"] <- slack_json[[message]]$latest_reply
    }
    #this will trigger for all child messages
    if (is.null(slack_json[[message]]$parent_user_id) == FALSE) { 
      messages_df[message, "ts_thread"] <- slack_json[[message]]$thread_ts
      messages_df[message, "parent_user_id"] <- slack_json[[message]]$parent_user_id
    }
  }
  
  return(messages_df)
}


#Run slack_json_to_dataframe() on all individual files in a channel (1 file / channel / day). 
#Bind them into a single df for each channel; add the channel name as a column.
#Finally, bind all of the individual channel dfs into a single dataframe for a given export!

#initialize the df for ALL THE MESSAGES across multiple days in all channels
all_channels_all_files_df <- setNames(data.frame(matrix(ncol = 12, nrow = 0)), 
                                      c("msg_id", "ts", "user", "type", "text",
                                        "reply_count", "reply_users_count", 
                                        "ts_latest_reply", "ts_thread, parent_user_id",
                                        "channel"))


for (channel in 1:length(channels_json)) {
  #initialize the df for ALL THE MESSAGES across multiple days in a single channel
    all_channel_files_df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), 
                             c("msg_id", "ts", "user", "type", "text",
                               "reply_count", "reply_users_count", 
                               "ts_latest_reply", "ts_thread, parent_user_id"))

  for (file_day in 1:length(channels_json[[channel]]$dayslist)) {
    #import the json file
      parentfolder_path <- paste0(slackexport_folder_path,"/",channels_json[[channel]]$name)
      filejson_path <- paste0(parentfolder_path, "/", channels_json[[channel]]$dayslist[[file_day]])
      import_file_json <-fromJSON(file = filejson_path)
    #initialize import_file_df for messages in a single day in a single channel
      import_file_df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)), 
                                 c("msg_id", "ts", "user", "type", "text",
                                   "reply_count", "reply_users_count", 
                                   "ts_latest_reply", "ts_thread, parent_user_id"))
    #convert json file to df
    import_file_df <- slack_json_to_dataframe(import_file_json)
    #bind the files together into a single df capturing all messages in a channel
    all_channel_files_df <- rbind(all_channel_files_df,import_file_df)
  }

  #Backfill channel name in the giant df for the all_channel_files_df you just created
  all_channel_files_df$channel <- channels_json[[channel]]$name
  #Bind all_channel_files together so all messages in all channels are in one file
  all_channels_all_files_df <- rbind(all_channels_all_files_df, all_channel_files_df)
}

#write the all files to a CSV in your R working directory
#format: exportfoldername_mindate_to_maxdate.csv
filename_mindate <- min(all_channels_all_files_df$ts) %>% as.numeric() %>% as.Date.POSIXct()
filename_maxdate <- max(all_channels_all_files_df$ts) %>% as.numeric() %>% as.Date.POSIXct()
  #Note exportfoldername was defined earlier before pulling in any of the files: exportname
slack_export_df_filename <- paste0(exportname,"_",filename_mindate,"_to_",filename_maxdate,".csv")
write.csv(all_channels_all_files_df, file = slack_export_df_filename)

#TODO - how does it handle orphaned threads? or deleted children? 
#TODO - make a users table with user metadata, write to csv
users_path <- paste0(slackexport_folder_path,"/users.json")
users_json <- fromJSON(file = users_path)
#initialize empty user df
user_list_df <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), 
                         c("user_id", "team_id", "name", "deleted", "real_name",
                           "tz", "tz_label", "tz_offset", "title", "display_name", 
                           "is_bot"))
#fill it with the appropriate fields from JSON
for (user in 1:length(users_json)) {
  user_list_df[user, "user_id"] <- users_json[[user]]$id
  user_list_df[user, "team_id"] <- users_json[[user]]$team_id
  user_list_df[user, "name"] <- users_json[[user]]$name
  user_list_df[user, "deleted"] <- users_json[[user]]$deleted
  #real_name is in a different place for bots - its nested in $profile
  if (is.null(users_json[[user]]$real_name) == FALSE) {
    user_list_df[user, "real_name"] <- users_json[[user]]$real_name
  }
  if (is.null(users_json[[user]]$profile$real_name) == FALSE) {
    user_list_df[user, "real_name"] <- users_json[[user]]$profile$real_name
  }
  user_list_df[user, "title"] <- users_json[[user]]$profile$title
  user_list_df[user, "display_name"] <- users_json[[user]]$profile$display_name
  user_list_df[user, "is_bot"] <- users_json[[user]]$is_bot
  #bots (?not sure who else) don't have time zone information. catch that null
  if (is.null(users_json[[user]]$tz) == FALSE) {
    user_list_df[user, "tz"] <- users_json[[user]]$tz
    user_list_df[user, "tz_label"] <- users_json[[user]]$tz_label
    user_list_df[user, "tz_offset"] <- users_json[[user]]$tz_offset
  }
  
}

#text data handling

#grab relevant libraries
library(dplyr) #pipe / data manipulation

#Read in the CSVs that were created in slack-export-data-handling.r
users_df <- user_list_df
channels_df <- channel_list
messages_df <- all_channels_all_files_df

#Add channel and user metadata into the all_files dataframe.
all_slack_data <- messages_df %>% left_join(users_df, by = c("user" = "user_id"))
all_slack_data <- all_slack_data %>% left_join(channels_df, by = c("channel" = "name"))

# Get questions
questions <- all_slack_data %>%
  filter(purpose == "Questions for class/worksheet on \"Dead or alive: Survival estimation\".") %>%
  filter(reply_count > 0) %>%
  select(text, ts, purpose) 
questions

# Get answers
answers <- all_slack_data %>%
  filter(purpose == "Questions for class/worksheet on \"Dead or alive: Survival estimation\".") %>%
  filter(is.na(reply_count)) %>%
  filter(!is.na(ts_thread)) %>%
  select(text, ts_thread)
answers

# Put everything together
all <- questions %>% 
  left_join(answers, by = c("ts" = "ts_thread")) %>%
  select(text.x, text.y) %>%
  rename("Question" = text.x,
         "Answer(s)" = text.y)

# Replace duplicates of questions
all[duplicated(all$Question),1] <- NA

# Make NA cells white so that we can't see them
all$Question <- cell_spec(all$Question, 
                          bold = TRUE,
                          color = ifelse(is.na(all$Question), "white", "black"))

# Display table
all %>%
  kbl(escape = FALSE) %>%
  kable_paper("hover", full_width = F)
