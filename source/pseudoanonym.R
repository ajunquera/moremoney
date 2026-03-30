# ...............................................................................
# ASSEGNO PER IL LAVORO - Pseudo-anonymization
# Author: Álvaro F. Junquera
# ...............................................................................

library(tidyverse)

# Reading data -------------
indi_ns_ss1 <- readRDS("intermediate/script01/indi_ns_ss1_190225.RDS")
indi_ns_ss2 <- readRDS("intermediate/script01/indi_ns_ss2_190225.RDS")

# Providers -----------
indi_ns_ente <- bind_rows(indi_ns_ss1, indi_ns_ss2) %>%
  select(axl_ente)

dict_ente <- data.frame(provider = unique(indi_ns_ente$axl_ente),
                        idprovider = 1:length(unique(indi_ns_ente$axl_ente)))

dict_ente$idprovider <- if_else(is.na(dict_ente$provider), NA, dict_ente$idprovider)

# Save -------------
saveRDS(dict_ente, "intermediate/dict_ente.RDS")
