
library("usethis")
usethis::use_git_config(user.name="villani", user.email= "giacomo.villani02@universitadipavia.it")
usethis::create_github_token() 
credentials::set_github_pat("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")
usethis::edit_r_environ()
usethis::git_sitrep()
