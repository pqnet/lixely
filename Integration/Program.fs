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
    let wf = dist {
        let! a = bernoulli 0.5
        let atime = if a then 1.0 else 3.0
        let! b = bernoulli 0.5
        let btime = if b then 1.0 else 3.0
        let! c = bernoulli 0.5
        let ctime = if c then 1.0 else 3.0

        let! d = bernoulli 0.5
        let dtime = (if d then 1.0 else 3.0) + (max atime btime)
        let! e = bernoulli 0.5
        let etime = (if e then 1.0 else 3.0) + btime
        let! f = bernoulli 0.5
        let ftime = (if f then 1.0 else 3.0) + (max ctime btime)
        return max (max dtime etime) ftime
    }
    //pdist wf

    let wf2 = dist {
        let! a = bernoulli 0.3
        let! b = bernoulli 0.7
        do! (a = b)
        return a
    }
    pdist wf2

//    let invokeRentalAgency =
//        dist {
//            let! success = bernoulli 0.85
//            if (success) then
//                let! cheap = bernoulli (0.7 / 0.85)
//                if (cheap) then
//                    return true,15.0
//                else
//                    return true,10.0
//            else
//                return false, 5.0
//        }
//    let invokeLoanAgent =
//        dist {
//            let! success = bernoulli 0.8
//            if (success) then
//                let! cheap = bernoulli (0.7 / 0.8)
//                if (cheap) then
//                    return true,15.0
//                else
//                    return true,10.0
//            else
//                return false, 5.0
//        }
//
//    let b = bernoulli 0.4
//    let c v = dist.Return v
//    let a = dist.Bind (b,c)
//    pdist a
    0
