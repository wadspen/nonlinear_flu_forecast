library(dplyr)

test_files <- list.files(pattern = "test_scores")
print(test_files)
all_test <- data.frame()

for (i in test_files) {
	file <- read.csv(i)
	print(head(file))
	all_test <- rbind(all_test, file)
     }

all_sum <- all_test %>%
	group_by(model, season, last_week, week_ahead) %>%
	summarise(mm = mean(mae), sdm = sd(mae), mc50 = mean(cover50),
		  sdc50 = sd(cover50), mc95 = mean(cover95), 
		  sdc95 = sd(cover95)) %>%
	mutate(sem = sdm/sqrt(1000), sec50 = sdc50/sqrt(1000),
	       sec95 = sdc95/sqrt(1000))


write.csv(all_sum, "samp_eval_sum.csv", row.names = FALSE)













