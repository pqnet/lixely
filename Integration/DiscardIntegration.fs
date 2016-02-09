module DiscardIntegration
open Swensen.Unquote.Extensions

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
    member dist.Spiattella(x:D<'T,_>):D<'T,_> = 
        let v,w = x {
            new FVectorSpace<_,Map<_,float>> with
            member self.Map v = Map.empty.Add (v,1.0)
            member self.Scale s x = Map.map (fun k v -> s * v) x
            member self.Add x y = x|> Map.fold (fun acc k v-> match acc.TryFind k  with None -> acc.Add(k,v)| Some w -> acc.Add(k,(v+w)) ) y
            member self.Zero = Map.empty
            }
        match w with
        | 0.0 -> fun v -> v.Zero,0.0
        | _ ->
            let l = v |> Map.toArray
            fun vspace ->
                let folder x (v,s) =
                    vspace.Map v
                    |>vspace.Scale s 
                    |>vspace.Add x
                Array.fold folder vspace.Zero l,w

    member inline dist.Run x = printf "." ; x// dist.Spiattella x
    member inline dist.Zero() : D<_,_> =
        fun fvspace -> (fvspace.Zero, 0.0)
    member inline dist.Return(x) : D<_,_>  =
        fun fvspace -> (fvspace.Map x, 1.0)

    member inline dist.Bind(x:D<_,_*float>, p: _ -> D<_,_>) : D<_,_> =
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

    member inline dist.cast (x:D<_,obj>): D<_,_> =
        fun fvspace ->
            let objvspace = {
                new FVectorSpace<_,obj> with
                member self.Map v = fvspace.Map v |> box
                member self.Scale s x = fvspace.Scale s (unbox x) |> box
                member self.Add x y = fvspace.Add (unbox x) (unbox y) |> box
                member self.Zero = fvspace.Zero |> box
            }
            let (r,w) = x objvspace
            unbox r, w
    member inline dist.ReturnFrom(x) =
        x

    [<CustomOperation("coerce", MaintainsVariableSpaceUsingBind=true)>]
    member inline dist.Coerce(x:D<_,_>, [<ProjectionParameter>] p) : D<_,_> =
        //let x = dist.Spiattella x
        dist { let! v = x in if p v then return v }

let dist = IntegratorBuilder

let bernoulli x : D<_,_> =
    fun fvspace ->
        (fvspace.Add (fvspace.Scale x (fvspace.Map true)) (fvspace.Scale (1.0-x) (fvspace.Map false)), 1.0)

let uniform ((z::l)) :D<_,_> =
    let l = List.toArray l
    let n = 1.0 / ((Array.length l |> float) + 1.0)
    fun fvspace ->
        let mutable r = fvspace.Map z
        for v in l do
            r <- fvspace.Add r (fvspace.Map v)
        fvspace.Scale n r,1.0

let UnionUniform<'T,'U> () :D<'T,'U> =
    let cases = Microsoft.FSharp.Reflection.FSharpType.GetUnionCases typeof<'T>
    cases |> Seq.map (fun x ->  Microsoft.FSharp.Reflection.FSharpValue.MakeUnion(x,[||])) |> Seq.cast<'T>|> Seq.toList |> uniform

let coin () = bernoulli 0.5

let d ()=
    dist {
        let! a = coin () //bernoulli 0.5
        let! b = coin () //bernoulli 0.5
        let! c = coin () // bernoulli 0.5
        coerce (a || b)
        coerce (c)
        return a
    }
let prob = d () floatVSpace
let prob2 = coin () floatVSpace
