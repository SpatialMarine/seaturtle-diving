# GAMM figs



library(ggplot2)
library(mgcv)
library(boot)
library(dplyr)
library(tidyr)

data<-read.csv("C:/Users/J. Menéndez Blázquez/SML_Dropbox/SML Dropbox/gitdata/seaturtle-diving/output/jess/plots/dive_summariesLabStates060224_FIN1304.csv")
data$datetime <- as.POSIXct(data$divetimedate,format="%Y-%m-%d %H:%M",tz="CET")
data$month<-format(data$datetime,"%m")

freq.df<-as.data.frame(table(data$organismID,data$month, data$states4))

freq.df <- freq.df %>% 
  rename("id"="Var1",
         "month"= "Var2",
         "state"="Var3")


prop_dives<-as.data.frame(freq.df %>%
                            group_by(month) %>% #id
                            mutate(prop =  Freq/sum(Freq)) %>% 
                            ungroup)

freq.df2<-as.data.frame(freq.df %>%
                          group_by(month,id) %>% #id
                          mutate(tot_divesMonth =sum(Freq)) %>% 
                          ungroup)

freq.df2<-as.data.frame(freq.df2 %>%
                          group_by(month,id) %>% 
                          mutate(dives_not = tot_divesMonth-Freq))


freq.df2<-as.data.frame(freq.df2 %>%
                          group_by(month,id) %>% #id
                          mutate(prop =  Freq/sum(Freq)) %>% 
                          ungroup)


freq.df2 <- freq.df2 %>% 
  rename("dives_performed"="Freq")

#dataframe now binomial

freq.df2$state<-as.factor(freq.df2$state)
freq.df2$id<-as.factor(freq.df2$id)
freq.df2$month<-as.numeric(freq.df2$month)


S1_data<-filter(freq.df2, state == "1")
S2_data<-filter(freq.df2, state == "2")
S3_data<-filter(freq.df2, state == "3")
S4_data<-filter(freq.df2, state == "4")
S1_data<- na.omit(S1_data)
S2_data<- na.omit(S2_data)
S3_data<- na.omit(S3_data)
S4_data<- na.omit(S4_data)


#ignore lag 0, only start form lag 1
#if sig. spikes waves in ACF chart, only 1/2 sig spikes use AR2


#check autocorrelation
library(tseries)
gamS1<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"),
           data=S1_data,
           method = 'REML',
           family = "binomial")
acf(resid(gamS1), main="acf(resid(m1))")
pacf(resid(gamS1), main="pacf(resid(m1))")

#lag 1
# state 2 dopes not need time stuff
gamS2<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"),
           data=S2_data,
           method = 'REML',
           family = "binomial")
acf(resid(gamS2), main="acf(resid(m1))")
pacf(resid(gamS2), main="pacf(resid(m1))")


gamS3<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"),
           data=S3_data,
           method = 'REML',
           family = "binomial")
acf(resid(gamS3), main="acf(resid(m2))")
pacf(resid(gamS3), main="pacf(resid(m2))")


gamS4<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"),
           data=S4_data,
           method = 'REML',
           family = "binomial")
acf(resid(gamS4), main="acf(resid(m3))")
pacf(resid(gamS4), main="pacf(resid(m3))")

plot(S4_data)

#using cc not bs as seasonality
##run gams for each state to assess seasonality

gamS1<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"), 
           data=S1_data,
           method = 'REML',
           family = "binomial")


crossing_df <- crossing(month = seq(1,12, length = 100),
                        id = levels(factor(S1_data$id)))



gamS1_preds_matrix <- predict(gamS1, 
                              newdata = crossing_df,
                              se.fit = TRUE,
                              re.form = NA, exclude=c("s(id)"))


gamS1_preds_df <- bind_cols(crossing_df,
                            gamS1_preds_matrix )


# ------------------------------------------------------------------------------
# Fig 6 ------------------------------------------------------------------------
# S1 plot -------------------------
p <- gamS1_preds_df %>% 
  ggplot(aes(x = month, y = inv.logit(fit)))+
  geom_jitter(data = S4_data,
              aes(y = prop),
              alpha = 0.15,
              size = 1.75,
              shape = 21,
              fill = "#EEAEEE",
              colour = "#CD96CD",
              width = 0.1,
              height = 0.01) +
  geom_ribbon(aes(ymin = inv.logit(fit - 1.96 * se.fit),
                  ymax = inv.logit(fit + 1.96 * se.fit)),
              alpha = 0.5,
              fill = "#ca7dcc" )+
  geom_line(color = "hotpink4") +
  labs(x = "Month", y="Proportion of time spent")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10, margin = margin(r = 8)),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(6, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

p

# export plot
p_png <- paste0(output_dir,"/fig/fig6a.png")
p_svg <- paste0(output_dir,"/fig/fig6a.svg")
ggsave(p_png, p, width=10, height=9, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=10, height=9, units="cm", dpi=300, bg="white")










# S2 --------------------------------------------------------------------------

gamS2<-gamm(cbind(dives_performed,dives_not) ~
              s(month, bs= "cc", k=6)+
              s(id, bs = "re"),correlation=corAR1(form=~1|month), 
            data=S2_data,
            method = 'REML',
            family = "binomial")

crossing_df <- crossing(month = seq(1,12, length = 100),
                        id = levels(factor(S2_data$id)))



gamS2_preds_matrix <- predict(gamS2$gam, 
                              newdata = crossing_df,
                              se.fit = TRUE,
                              re.form = NA, exclude=c("s(id)"))


gamS2_preds_df <- bind_cols(crossing_df,
                            gamS2_preds_matrix )

# plot ---------------------------

p <-gamS2_preds_df %>% 
  ggplot(aes(x = month, y = inv.logit(fit)))+
  geom_jitter(data = S4_data,
              aes(y = prop),
              alpha = 0.25,
              size = 1.75,
              shape = 21,
              fill = "khaki2",
              colour = "#CDC673",
              width = 0.1,
              height = 0.01) +
  geom_ribbon(aes(ymin = inv.logit(fit - 1.96 * se.fit),
                  ymax = inv.logit(fit + 1.96 * se.fit)),
              alpha = 0.55,
              fill = "#E69F00" )+
  geom_line(color = "#8B5A00") +
  labs(x = "Month", y="Proportion of time spent")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10, margin = margin(r = 8)),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(6, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")
p

# export plot
p_png <- paste0(output_dir,"/fig/fig6b.png")
p_svg <- paste0(output_dir,"/fig/fig6b.svg")
ggsave(p_png, p, width=10, height=9, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=10, height=9, units="cm", dpi=300, bg="white")





# S3 Plot ----------------------------------------------------------------------

gamS3<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"), 
           data=S3_data,
           method = 'REML',
           family = "binomial")




crossing_df <- crossing(month = seq(1,12, length = 100),
                        id = levels(factor(S3_data$id)))



gamS3_preds_matrix <- predict(gamS3, 
                              newdata = crossing_df,
                              se.fit = TRUE,
                              re.form = NA, exclude=c("s(id)"))


gamS3_preds_df <- bind_cols(crossing_df,
                            gamS3_preds_matrix )


p <-gamS3_preds_df %>% 
  ggplot(aes(x = month, y = inv.logit(fit)))+
  geom_jitter(data = S4_data,
              aes(y = prop),
              alpha = 0.15,
              size = 1.75,
              shape = 21,
              fill = "skyblue1",
              colour = "skyblue2",
              width = 0.1,
              height = 0.01) +
  geom_ribbon(aes(ymin = inv.logit(fit - 1.96 * se.fit),
                  ymax = inv.logit(fit + 1.96 * se.fit)),
              alpha = 0.55,
              fill = "skyblue3" )+
  geom_line(color = "#104E8B") +
  labs(x = "Month", y="Proportion of time spent")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10, margin = margin(r = 8)),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(6, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

p

# export plot
p_png <- paste0(output_dir,"/fig/fig6c.png")
p_svg <- paste0(output_dir,"/fig/fig6c.svg")
ggsave(p_png, p, width=10, height=9, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=10, height=9, units="cm", dpi=300, bg="white")








# S4 plot ----------------------------------------------------------------------

gamS4<-gam(cbind(dives_performed,dives_not) ~
             s(month, bs= "cc", k=6)+
             s(id, bs = "re"), 
           data=S4_data,
           method = 'REML',
           family = "binomial")

crossing_df <- crossing(month = seq(1,12, length = 100),
                        id = levels(factor(S4_data$id)))

gamS4_preds_matrix <- predict(gamS4, 
                              newdata = crossing_df,
                              se.fit = TRUE,
                              re.form = NA, exclude=c("s(id)"))

gamS4_preds_df <- bind_cols(crossing_df,
                            gamS4_preds_matrix )

# State 4 plot
p <-gamS4_preds_df %>% 
  ggplot(aes(x = month, y = inv.logit(fit)))+
  geom_jitter(data = S4_data,
              aes(y = prop),
              alpha = 0.15,
              size = 1.75,
              shape = 21,
              fill = "coral",
              colour = "coral2",
              width = 0.1,
              height = 0.01) +
  geom_ribbon(aes(ymin = inv.logit(fit - 1.96 * se.fit),
                  ymax = inv.logit(fit + 1.96 * se.fit)),
              alpha = 0.65,
              fill = "coral" )+
  geom_line(color = "darkred") +
  labs(x = "Month", y="Proportion of time spent")+
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))+
  # theme
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 10, margin = margin(r = 8)),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks = element_line(size = 0.7),
        axis.ticks.length = unit(6, "pt"),
        panel.border = element_rect(color = "black", fill = NA, size = 1.2),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none")

p

# export plot
p_png <- paste0(output_dir,"/fig/fig6d.png")
p_svg <- paste0(output_dir,"/fig/fig6d.svg")
ggsave(p_png, p, width=10, height=9, units="cm", dpi=300, bg="white")
ggsave(p_svg, p, width=10, height=9, units="cm", dpi=300, bg="white")


  







s1_graph
s2_graph
s3_graph
s4_graph


summary(gamS1)
summary(gamS2$gam)
summary(gamS3)
summary(gamS4)

