
//distribuzione: una funzione che prende una tupla di valori casuali, e restituisce un valore.

//sampler: un oggetto che tiene traccia di tutte le variabili aleatorie create, e di quanta entropia e' necessaria per fare un sample di tutte.

//variabile aleatoria: un oggetto con una distribuzione, ed eventualmente correlato a variabili aleatorie precedenti.
module Distributions

//value type for generator (and thus for uniform) samples
type Gen_Sample_T = double

//abstract type for random source generator:
type IGenerator =
    abstract NextValue: Gen_Sample_T*IGenerator

// functional monadic-style wrapper for generator 
let rec genRaw status =
    let randomFunc status =
        let (g:RandomTools.RandomSource),idx = status
        let newStatus = (g,idx+1)
        let value = g.GetSample(idx)
        (value,newStatus)
    {
        new IGenerator with
            member g.NextValue = let value,status = randomFunc status
                                 (value,genRaw status)
    }

// new generator constructor
let gen () =
    genRaw ((new RandomTools.RandomSource () ) ,0)

let fastGen =
    let random = new System.Random()
    {
        new IGenerator with
            member g.NextValue = random.NextDouble(),g
    }
type IDistribution<'T> =
    abstract NextValue: IGenerator -> 'T * IGenerator

(*type IDiscreteDistribution<'T> =
    inherit IDistribution<'T>
    abstract Domain: 'T seq
    abstract GetProb: 'T -> float
*)
type DiscreteDistribution<'T>(domain:'T array, probabilities:float array) =
    do assert (domain.Length = probabilities.Length)
    let c = Seq.scan (+) 0.0 probabilities|> Seq.skip 1 |> Seq.toArray
    let normalize =  Array.map (fun x -> x/ c.[c.Length - 1])
    let cumul = normalize c
    let probs = normalize probabilities 
    member d.probabilities = probs
    member d.domain = domain
    member d.cumulative = cumul
    interface IDistribution<'T> with
        member d.NextValue g =
            let v,g1 = g.NextValue
            let idx =abs(array.BinarySearch(cumul,v) + 1)
            domain.[idx],g1
type ConstantDistribution<'T> (value:'T) =
    inherit DiscreteDistribution<'T>([|value|],[|1.0|])
    interface IDistribution<'T> with
        member d.NextValue g = value,g
    member d.Value = value


type DistributionBuilder () =
        
    //member d.Bind ((expr:IDiscreteDistribution<'T>), (f: 'T -> IDistribution<'U>))
    member d.Bind ((expr:ConstantDistribution<'T>), (f: 'T -> 'V when 'V :> IDistribution<'U> )) =
        f expr.Value
    member d.Bind ((expr:DiscreteDistribution<'T>), (f: 'T -> 'V when 'V :> DiscreteDistribution<'U>)) =
        let domain,probab = 
            expr.domain
            |> Seq.map f
            |> Seq.map (fun x -> Seq.zip x.domain x.probabilities)
            |> Seq.concat
            |> Seq.groupBy fst
            |> Seq.map (fun (k,s) -> k,Seq.map snd s |> Seq.fold (+) 0.0 )
            |> Seq.toArray
            |> Array.unzip
        DiscreteDistribution(domain,probab)


    member d.Bind ((expr:IDistribution<'T>), (f: 'T -> 'V when 'V :> IDistribution<'U>)) =
        {
            new IDistribution<'U> with
                member d.NextValue g =
                    let v,g1 = expr.NextValue g
                    (f v).NextValue g1
        }
    member d.Return (expr:'T) =
        ConstantDistribution(expr)

    member d.ReturnFrom (expr:'V when 'V :> IDistribution<'T>) =
        expr


let dist = new DistributionBuilder()

module Dist =
    let rec getSampleSeq (d:IDistribution<'T>) (g:IGenerator) =
        Seq.unfold (fun gen -> Some (d.NextValue gen)) g

    let uniform = {
        new IDistribution<double> with
                    member d.NextValue g = g.NextValue
        }
    let bernoulli p = DiscreteDistribution([|true;false|],[|p;1.0-p|])

    let normal m s =
        let gaussianBoxMuller (m:float<_>) (sigma:float<_>) =
            dist {
                let! u = uniform
                let! v = uniform
                return m + sigma * sqrt(-2. * log(1.-u)) * cos(2. * System.Math.PI * v)
            }
        gaussianBoxMuller m s