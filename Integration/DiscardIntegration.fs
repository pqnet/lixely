module Integration

let private pairwise (f1,f2) (a1:'T,a2:'U) :'V*'W= (f1 a1,f2 a2)

type D<'T> =
    //this is private
    | Integrator of (('T -> obj) -> (float->obj->obj) -> (obj -> obj -> obj) -> obj -> obj option)
    //              function           scaling            aggregation          zero   result
    member x.intfun<'U> (smap:'T->'U) (smul:float->'U->'U) (ssum:'U->'U->'U) (zero:'U) =
        let (Integrator(f)) = x
        let omap = smap >> box
        let omul = fun s -> unbox >> (smul s) >> box
        let osum = fun a-> unbox >> (unbox a |> ssum) >> box
        let ozero = box zero
        f omap omul osum ozero |> Option.map unbox<'U>
    static member Naught = Integrator (fun f m a z -> None) : D<'T>
        
let inline private ret (x:'T) =
    let q f m a z =
        f x |> Some
    q |> Integrator

let inline private fmap (x: 'T ->'U) (d:D<'T>) =
    let if1 f (*s a z*) = d.intfun (x>>f) (* s a z*)
    Integrator if1

let inline private join (d:D<D<'T>>)(*: D<'T> *)=
    let q f m a z =
        //outf :: D<'T> -> 'T integrando f
        let outf (d:D<'T>)  = d.intfun f m a z |> Option.map (fun x -> x,1.0)
        let outm s = Option.map (fun (x,wx) -> m s x,wx * s)
        let outa x y =
            match (x,y) with
            | (Some (x,wx), Some (y,wy)) -> Some (a x y,(wx+wy))
            | (None,None) -> None
            | (Some _, None) -> x 
            | (None,Some _) -> y
        let outz = None
        let v = d.intfun outf outm outa outz
        match v with
        | None -> None
        | Some None -> None
        | Some (Some (x,w)) -> m (1.0/w) x |> Some
    q |> Integrator

type IntegratorBuilder =
    | IntegratorBuilder
    member inline x.Zero<'T> () = D.Naught :D<'T>
    member inline x.Return(v:'T) = ret v
    member inline x.Bind(v:D<'T>,c:'T -> D<'U>) =
        fmap c v |> join

    [<CustomOperation("coerce",MaintainsVariableSpaceUsingBind=true)>]
    member inline x.Coerce(p1:D<'T>,[<ProjectionParameter>] p2:'T->bool):D<'T> =
        x.Bind (p1,fun e -> if p2 e then x.Return e else x.Zero ())

let dist = IntegratorBuilder
let bernoulli x = Integrator <| fun f s a z -> a (s x (f true)) (s (1.0-x) (f false)) |> Some

let b = true || false
let d = 
    dist {
        let! a = bernoulli 0.5
        let! b = bernoulli 0.5
        coerce (a)
        return a
    }

let prob = d.intfun (fun x -> if x then 1.0 else 0.0) (*) (+) 0.0
