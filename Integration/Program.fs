// Ulteriori informazioni su F# all'indirizzo http://fsharp.net
// Per ulteriori informazioni, vedere il progetto 'Esercitazione su F#'
module Main
open Integration

let pdist (d:D<'T>) =
    let mapper (x:'T) = Seq.singleton (1.0,x)
    let multp s x = Seq.map (fun (x,(y:'T)) -> (x*s,y)) x
    let sump = Seq.map fst >> Seq.reduce (+)
    let collect a b =
        Seq.append a b 
        |> Seq.groupBy snd 
        |> Seq.map (fun (k,v) -> sump v,k) 
    d.intfun 
        <| mapper
        <| multp
        <| collect
    |> Seq.iter (printf "%A") 


[<EntryPoint>]
let main argv =
    let b = bernoulli 0.4
    let c v = dist.Return v
    let a = dist.Bind (b,c)
    pdist a
    0
