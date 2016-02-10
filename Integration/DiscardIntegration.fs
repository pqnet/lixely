module DiscardIntegration
//open Swensen.Unquote.Extensions

let inline ( *=| ) (x: _ byref) s = (^U : (member Scale: float -> unit) (x, s))
let inline ( +=| ) (x: _ byref) y = (^U : (member Add: 'U -> unit) (x, y))

let inline origin () = (^U : (static member Origin: unit -> ^U) ())

let inline make<'T when ^T : (new : unit -> ^T)>() = new 'T()

[<StructuredFormatDisplay("F{V}")>]
type FloatVector() =
  [<DefaultValue>] val mutable v : float

  static member inline Origin () =
      let mutable r = make<FloatVector>()
      r.v <- 0.0
      r

  member inline x.V = x.v
  member inline x.Scale (s:float) = x.v <- x.v * s
  member inline x.Add (y:FloatVector) = x.v <- x.v + y.v

[<StructuredFormatDisplay("[{Map}]")>]
type MapVector<'T when 'T : comparison>() =
  [<DefaultValue(false)>] val mutable map : Map<'T, float>

  static member Origin () =
      let mutable r = make<MapVector<'T>>()
      r.map <- Map.empty
      r

  member x.Map = x.map
  member x.Scale s = x.map <- Map.map (fun k v -> v * s) x.Map
  member x.Add (y:MapVector<'T>) =
      let combine map k v =
          match Map.tryFind k map with
          | None -> map.Add(k, v)
          | Some w -> map.Add(k, v+w)
      x.map <- Map.fold combine x.Map y.Map


[<StructuredFormatDisplay("{Stringed}")>]
type DistResult<'U when ^U : (member  Scale : float -> unit)
                    and ^U : (member    Add :    ^U -> unit)
                    and ^U : (static member Origin :  unit -> ^U)>() =
  [<DefaultValue(false)>] val mutable v : 'U
  [<DefaultValue>] val mutable w : float

  member inline x.Stringed = sprintf "[%A : %A]" x.v x.w

  static member inline Origin () =
      let mutable r = make<DistResult<_>>()
      r.v <- origin ()
      r.w <- 0.0
      r

  member inline x.Scale s =
      &x.v *=| s
      x.w <- x.w * s

  member inline x.Add (y:DistResult<'U>) =
      &x.v +=| y.v
      x.w <- x.w + y.w


type D<'T,'U  when ^U : (       member  Scale : float -> unit)
               and ^U : (       member    Add :    ^U -> unit)
               and ^U : (static member Origin :  unit -> ^U  )> = ('T -> 'U) -> DistResult<'U>


let inline bool2float b =
    let mutable r = make<FloatVector>()
    r.v <- if b then 1.0 else 0.0
    r

let inline zero () = origin ()

let inline result v =
    let mutable r = make<DistResult<_>>()
    r.v <- v
    r.w <- 1.0
    r

let inline makeMap v =
    let mutable r = MapVector<_>.Origin()
    r.map <- r.map.Add(v, 1.0)
    r

type IntegratorBuilder =
    | IntegratorBuilder
    
    member inline dist.Spiattella(x:D<'T,MapVector<_>>) : D<'T,_> =
        let r = x makeMap
        match r.w with
        | 0.0 -> fun f -> zero ()
        | w ->
            let a = Map.toArray r.v.map
            fun f ->
                let mutable r : DistResult<_> = zero ()
                r.w <- w
                for (v,s) in a do
                    let mutable tmp = f v
                    &tmp *=| s
                    &r.v +=| tmp
                r

    member inline dist.Run x = x // dist.Spiattella x
    member inline dist.Zero() : D<_,_> = fun f -> zero ()
    member inline dist.Return(x) : D<_,_>  = fun f -> result (f x)
    member inline dist.ReturnFrom(x) : D<_,_> = x

    member inline dist.Bind(x:D<_,_>, p: _ -> D<_,_>) : D<_,_> =
        fun f ->
            let outf v = p v f
            let mutable r = x outf
            match r.w with
            | 0.0 -> zero ()
            | w -> &r.v *=| (1.0 / w) ; r.v

    [<CustomOperation("coerce", MaintainsVariableSpaceUsingBind=true)>]
    member inline dist.Coerce(x, [<ProjectionParameter>] p) =
        dist { let! v = x in if p v then return v }

let dist = IntegratorBuilder

let inline bernoulli x f =
        let mutable ftrue = f true
        let mutable ffalse = f false
        &ftrue *=| x
        &ffalse *=| (1.0 - x)
        &ftrue +=| ffalse
        result ftrue

let inline uniform l : D<_,_> =
    match l with
    | [] -> dist.Zero()
    | (hd::tl) ->
        let a = List.toArray tl
        let s = 1.0 / (1 + Array.length a |> float)
        fun f ->
            let mutable r = f hd
            for v in a do
                &r +=| f v
            &r *=| s                
            result r

let inline UnionUniform () : D<'T,_> =
    Microsoft.FSharp.Reflection.FSharpType.GetUnionCases typeof<'T>
    |> Array.map (fun x -> Microsoft.FSharp.Reflection.FSharpValue.MakeUnion(x, [||]) :?> _)
    |> Array.toList
    |> uniform

let inline coin () = bernoulli 0.5

let inline myprob () =
   dist {
        let! a = coin ()
        let! b = coin ()
        let! c = coin ()
        coerce (a || b)
        coerce (c)
        return a
    } <| bool2float

let eval_bernoulli x = bernoulli x bool2float

printfn "%A" (eval_bernoulli 0.7)
printfn "%A" (myprob ())

let mutable v = FloatVector.Origin()
printfn "%A" v
v.v <- 12.0
printfn "%A" v
&v +=| v
printfn "%A" v
&v *=| 3.2
printfn "%A" v
