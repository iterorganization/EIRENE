MODULE json_module
    type,public :: json_value
    end type
    type,public :: json_core
    end type
      
    integer,parameter,public :: json_RK = 0
    integer,parameter,public :: json_IK = 0
    integer,parameter,public :: json_LK = 0
    integer,parameter,public :: json_CK = selected_char_kind("DEFAULT")
    integer,parameter,public :: json_CDK = selected_char_kind("DEFAULT")

END MODULE
