module Integration


type D<'T> =
    | Integrator of (('T -> obj) -> (float->obj->obj) -> (obj -> obj -> obj) -> obj)
    //                function           scaling            aggregation        result
    member x.intfun<'U> (smap:'T->'U) (smul:float->'U->'U) (ssum:'U->'U->'U) =
        let (Integrator(f):D<'T>) = x
        let omul = fun s -> unbox >> (smul s) >> box
        let osum = fun a-> unbox >> (unbox a |> ssum) >> box
        f (smap>>box) omul osum |> unbox<'U>
        
let inline dni (x:'T) = Integrator <| fun f s a -> s 1.0 (f x)

let inline D (x: 'T ->'U) (d:D<'T>) =
    let if1 f s a = d.intfun (x>>f) s a
    Integrator if1

//troiaio. Ma che stavo facendo?
let dmi2 (Integrator(if1):D<D<'T>>): D<'T> =
    let sumD (Integrator(if1):D<'T>) (Integrator(if2):D<'T>) : D<'T> =
        fun f s a ->
            let v1 = if1 f s a
            let v2 = if2 f s a
            a v1 v2
        |> Integrator
    let scaleD (x:float) (Integrator(if1):D<'T>) : D<'T> =
        fun f s a ->
            s x <| if1 f s a
        |> Integrator
    let sumB (o1) (o2) =
        sumD (unbox o1:D<'T>) (unbox(o2)) |> box
    let scaleB (sv) (o) =
        scaleD sv (unbox(o):D<'T>) |> box
    fun f s a -> if1 (fun (Integrator(g)) -> g f s a ) (scaleB) (sumB)
    |> Integrator

let inline dmi (d:D<D<'T>>): D<'T> =
    fun f s a ->
        d.intfun (fun d -> d.intfun f s a) s a
    |> Integrator

//non esiste la possibilita' che qualcosa sia moltiplicabile per scalare. Ma questo non rompe le unita' di misura?
(*
let inline dehlay (Integrator(f): D< ^T > ) () =
    let boxsum (a:obj) (b:obj) = (unbox a: ^T) + (unbox b : ^T) |> box
    let boxmult (a:float) (b:obj) = (a:float) * (unbox b : ^T) |> box
    f box boxmult boxsum |> unbox< ^T > 
*)

type IntegratorBuilder =
    | IntegratorBuilder
    member inline x.Return(v:'T) = dni v
    member inline x.Bind(v:D<'T>,c:'T -> D<'U>) =
        let a = D c v 
        dmi a
    member inline x.Bind(v:bool,c:unit->D<'T>) =
        let cg = if v then x.Bind(c (),fun v -> x.Return(Some v)) else x.Return None
        Integrator <| fun f s a ->
            let ff = Option.map (fun x -> 1.0,f x)
            let ss (scale:float) :(float*obj) option -> _ = Option.map <| fun (p,v) -> p * scale,s scale v
            let aa p1 p2 =
                match p1,p2 with
                | None,x -> x
                | Some x,None -> Some x
                | Some (p1:float,v1), Some (p2,v2) -> Some (p1+p2,a v1 v2)
            match cg.intfun ff ss aa with
            | Some (p,v) -> s (1.0/p) v
            | None -> failwith "the condition imposed is almost never satisfied for the given priors"
//            let cont = x.Bind(v,fun guard -> if guard then x.Bind(c(),Some()) else None)
//            let support = ( v.intfun (fun b -> if b then 1.0 else 0.0) (*) (+))
//            v.intfun (fun b -> cont.intfun (fun x -> s (if b then 1.0 else 0.0) <| f x) s a) 
   
let dist = IntegratorBuilder
let bernoulli x = Integrator <| fun f s a -> a (s x (f true)) (s (1.0-x) (f false))
