module Integration

type FVectorSpace<'T,'U> =
    abstract member Map : 'T -> 'U
    abstract member Scale : float -> 'U -> 'U
    abstract member Add : 'U -> 'U -> 'U
    abstract member Zero : 'U

let floatVSpace = {
        new FVectorSpace<bool,float> with
        member this.Map v = if v then 1.0 else 0.0
        member this.Scale s x = s * x
        member this.Add a b = a + b
        member this.Zero = 0.0
    }

type D<'T,'U> = FVectorSpace<'T,'U> -> 'U*float

type IntegratorBuilder =
    | IntegratorBuilder
    member inline this.Zero() : D<_,_> =
        fun fvspace -> (fvspace.Zero, 0.0)
    member inline this.Return(x) : D<_,_>  =
        fun fvspace -> (fvspace.Map x, 1.0)
    member inline this.Bind(x:D<_,_>, p: _ -> D<_,_>) : D<_,_> =
        fun fvspace ->
        let outfvs =
            {
                new FVectorSpace<_,_> with
                member self.Map v = p v fvspace
                member self.Scale s ((x,w)) = fvspace.Scale s x,s * w
                member self.Add ((x,wx)) ((y,wy)) = fvspace.Add x y,(wx+wy)
                member self.Zero = fvspace.Zero, 0.0
            }
        match x outfvs with
        | (_, 0.0) -> outfvs.Zero
        | (r, w)   -> outfvs.Scale (1.0/w) r

    [<CustomOperation("coerce", MaintainsVariableSpaceUsingBind=true)>]
    member inline this.Coerce(p1, [<ProjectionParameter>] p2) : D<_,_> =
        IntegratorBuilder { let! v = p1 in if p2 v then return v }

let dist = IntegratorBuilder

let bernoulli x : D<_,_> =
    fun fvspace ->
        (fvspace.Add (fvspace.Scale x (fvspace.Map true)) (fvspace.Scale (1.0-x) (fvspace.Map false)), 1.0)

let coin () = bernoulli 0.5

let d =
    dist {
        let! a = coin () //bernoulli 0.5
        let! b = coin () //bernoulli 0.5
        let! c = coin ()// bernoulli 0.5
        coerce (a || b)
        coerce (c)
        return a
    }



let prob = d floatVSpace
let prob2 = coin () floatVSpace

//[<EntryPoint>]
let main args =
  printfn "%A" prob
  0
