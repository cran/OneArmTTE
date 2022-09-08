# These global variables are declared to eliminate associated R cmd check warnings.
# There is no other identified functional impact of these global declarations.

utils::globalVariables(
  c('.',
    'censor',
    'cte',
    'dropoutRate',
    'dropoutTime',
    'enrollTime',
    'fail',
    'failTime',
    'duration',
    'enrollTime',
    'finish',
    'lambda',
    'N',
    'naive_flag',
    'origin',
    'period',
    'rate',
    'time',
    't.cut'
  )
)
