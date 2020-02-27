# ? Use a macro @PowerNode to generate it ?

abstract type Node()
    """
    # ! A node is defined by a:
    #       - label (can be a name or a number)
    #       - type (Bus or Component)
    #       - Input flows
    #       - Output flows
    """

    subtype Bus()
    """
        # ? define as function, as (mutable) struct or as an abstract type?
        # ? What are the main characteristics for a bus
        # ? Use Parameters.jl to define these characteristics
    """
    end

    subtype Component()
    """
        # ? define as function, as (mutable) struct or as an abstract type?
        # ? Use Parameters.jl to define these characteristics
        # ? Storage or no storage

        # ! Must have one or more outputs
        # ! Each component has general characteristics, as defined in the Power Node Framework model
        #
    """
    end

end
