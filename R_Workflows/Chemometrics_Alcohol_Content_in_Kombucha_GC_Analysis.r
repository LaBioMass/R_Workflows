# Upload packages 
library(ggplot2)
library(dplyr)
library(car)  
library(MASS) 

# Create a data frame
dados <- data.frame(
  Temperatura = c(110, 110, 110, 110, 110, 110, 70, 110, 70, 150, 150, 110, 150, 110, 110, 70, 53.43, 53.43, 110, 150, 110, 110, 166.57, 110, 70),
  Tempo = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 5.0, 3.5, 2.0, 2.0, 2.0, 5.62, 5.0, 3.5, 1.38, 5.0, 3.5, 3.5, 1.38, 5.0, 5.62, 3.5, 3.5, 3.5, 2.0),
  Ultrassom = factor(c("sem", "sem", "com", "com", "com", "com", "sem", "sem", "com", "sem", "com", "sem", "com", "sem", "com", "com", "com", "sem", "sem", "sem", "com", "sem", "com", "com", "sem")),
  Intensidade = c(1.02897E+07, 1.07284E+07, 1.10983E+07, 1.10926E+07, 1.16278E+07, 1.15247E+07, 1.0985E+06, 1.00828E+07, 1.06048E+06, 1.25084E+06, 9.97403E+06, 1.01945E+07, 2.9545E+06, 1.01836E+07, 8.22732E+06, 1.33499E+06, 5.86538E+06, 1.00367E+06, 1.09216E+07, 89765, 1.15243E+06, 1.05246E+07, 363319, 1.11052E+07, 1.13774E+06))

#You can also upload direct from your directory

# Models (to creat a model the fuction is lm)

modelo_linear <- lm(Intensidade ~ Temperatura + Tempo + Ultrassom, data = dados)
modelo_quadratico <- lm(Intensidade ~ Temperatura + I(Temperatura^2) + Tempo + I(Tempo^2) + Ultrassom, data = dados)
modelo_cubico <- lm(Intensidade ~ Temperatura + I(Temperatura^2) + I(Temperatura^3) + 
                      Tempo + I(Tempo^2) + I(Tempo^3) + Ultrassom, data = dados)

# Models (resume)
summary(modelo_linear)
summary(modelo_quadratico)
summary(modelo_cubico)

# In this case, we use the quadratic model 

# Coeficients of quadratic models
beta_1 <- coef(modelo_quadratico)["Temperatura"]
beta_2 <- coef(modelo_quadratico)["I(Temperatura^2)"]
beta_3 <- coef(modelo_quadratico)["Tempo"]
beta_4 <- coef(modelo_quadratico)["I(Tempo^2)"]

# Derivatives for temperature and time
T_otimo <- -beta_1 / (2 * beta_2)
P_otimo <- -beta_3 / (2 * beta_4)

# Optimized values
T_otimo
P_otimo

#Length of values (determined by your data)
temp_seq <- seq(from = 50, to = 150, length.out = 20)  
tempo_seq <- seq(from = 1, to = 5, length.out = 20)  

# Create a grid of Temperature and Time combinations
grade <- expand.grid(Temperatura = temp_seq, Tempo = tempo_seq)

# Exclude any compound that is insignificant
grade$Ultrassom <- factor("sem", levels = c("sem", "com"))

# Response for significant terms
grade$Intensidade <- predict(modelo_quadratico, newdata = grade)

# Surface graph
library(ggplot2)
library(viridis)

ggplot(grade, aes(x = Temperatura, y = Tempo, z = Intensidade)) +
  geom_tile(aes(fill = Intensidade), alpha = 12.5) +  
  stat_contour(aes(z = Intensidade), color = "black") +  
  scale_fill_viridis_c() +  
  theme_minimal() +
  labs(title = "Surface Response",
       x = "Temperature (°C)",
       y = "Time (h)",
       fill = "Intensidade") +
geom_point(aes(x = T_otimo, y = P_otimo), color = "red", size = 3) + 
  annotate("text", x = T_otimo, y = P_otimo, label = '', color = "red", vjust = -1)

#Interaction between significant terms
ggplot(grade, aes(x = Temperatura, y = Intensidade, color = factor(Tempo))) +
  geom_line(size = 1.0) +
  labs(title = "Interaction between temperature and time", x = "Temperature (°C)", y = "Intensity", color = "Time") +
  theme_bw()
