parameterGroup <- function(a_list, b_list, c_list){
    "
    Example:
        res = parameterGroup(c(1,2,3), c(1,2,3), c(1,2,3))
    "
    # list length
    a_len = length(a_list)
    b_len = length(b_list)
    c_len = length(c_list)
    # group parameters
    par_groups = list()
    cnt = 1
    for (a in a_list) {
       for (b in b_list) {
          for(c in c_list){
              par_groups[[cnt]] = c(a,b,c)
              cnt = cnt + 1
          }
       }
    }
    return (par_groups)
}