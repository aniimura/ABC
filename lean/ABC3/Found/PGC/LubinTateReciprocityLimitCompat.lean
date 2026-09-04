import ABC3.Found.PGC.LubinTateGeneratorSequence

/-!
# `reciprocityMap` を無限compatible列の上で評価した値のn跨ぎ両立性(`sorry` 無し)

節目(5)(射影極限 `Gal(L_π/K)≅𝒪_K^×`)の最終組み立ての核心部分: 前段
(`Found/PGC/LubinTateTowerCompatible.lean::reciprocityMap_pred_eq_
map_succ`)と`Found/PGC/LubinTateGeneratorSequence.lean::psiGenSeq`
(無限compatible列)を実際に繋ぎ合わせ、「同じ`σ`について、`reciprocity
Map`を無限列`(psiGenSeq n).pt`の上で評価した値が射影系をなす」ことを
確立する。

## 鍵となる補題

- `reciprocityMap_congr`: `reciprocityMap`は評価点`x`が(値として)
  同じであれば、どの証明項を使っても同じ値を返す——`x=x'`から
  `subst`一発で`rfl`により従う、`Prop`の証明無関係性からの直接の
  帰結。
- `reciprocityMapLimitCompat`: 上の2つを組み合わせて実際にn跨ぎの
  両立性を示す。`psiGenStepResult`の`hcompat`(既出、`π·(次の生成元)
  =(前の生成元)`)を使い、`reciprocityMap_pred_eq_map_succ`の要求する
  `y`(`hyψ`/`hyn`/`hymem`)を`rw[g.hcompat]`で`(psiGenSeq m)`の対応
  する事実へ変換してから適用し、最後に`reciprocityMap_congr`で
  `psiGenStepResult`基準の表現から`psiGenSeq`基準の表現へ載せ替える。

## 踏んだ罠(`tools/lean-idioms.md` #34 の続報)

`reciprocityMapLimitCompat`の返す型を**省略**した版(`noncomputable
def`)は、`Eq.trans`という単純な操作にもかかわらず「型を独立に書いて
いない」ことが災いし、`?m`という未解決のメタ変数を残したまま`exact`
が失敗する(こちらは`Classical.choice`由来の罠とは**別の**、単に
「型注釈なしの`def`+`by`ブロック」が大きすぎる証明で型推論に失敗する
という、より平凡な現象)。★対処: 今回は**型注釈を省略せず**、むしろ
`FiniteDimensional`インスタンスを2つ明示的な`[...]`引数として追加
する(`reciprocityMap`が要求するインスタンスが`.hfd`経由でしか
手に入らないため)ことで、独立に書いた型注釈が今度は問題なく通った
——`Classical.choice`由来の罠(#34)は「型を独立に書くと落ちる」
だったが、今回は逆に「型を書かないと(この特定の証明の大きさでは)
推論できない」という、状況に応じた使い分けが要ることを実測した。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-- **`reciprocityMap`は評価点の値だけで決まる**——証明項(`hxψなど`)
が違っても、評価点`x`自身が(値として)一致すれば同じ値を返す。
`Prop`の証明無関係性から`subst`+`rfl`で直ちに従う。 -/
theorem reciprocityMap_congr {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x x' : K.closure) (h : x = x')
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (hxψ' : x' ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn' : x' ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem' : x' ∈ IntermediateField.adjoin K.carrier ({x'} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x'} : Set K.closure))]
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ =
      reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x' hxψ' hxn' hmem' σ := by
  subst h
  rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`reciprocityMap`を無限compatible列の上で評価した値のn跨ぎ両立性**:
同じ大域的`σ`について、`(psiGenSeq(m+1)相当)`で評価した`reciprocity
Map`を`principalUnits K π (m+1)`まで落としたものが、`(psiGenSeq m)`
で評価した`reciprocityMap`に一致する。`reciprocityMap_pred_eq_map_
succ`(既出)を`psiGenStepResult`が与える具体的な生成元・compat性
(`hcompat`)へ適用し、`reciprocityMap_congr`で表現を`psiGenSeq`基準
へ載せ替えるだけ。 -/
theorem reciprocityMapLimitCompat {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (m : ℕ)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt} :
          Set K.closure))]
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure))]
    (σ : K.closure ≃ₐ[K.carrier] K.closure) :
    (QuotientGroup.map (principalUnits K π (m + 1 + 1)) (principalUnits K π (m + 1)) (MonoidHom.id _)
        (principalUnits_succ_le K π (m + 1)))
      (reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) (by omega)
        (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).pt
        (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hψ
        (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hn
        (psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m)).hmem σ) =
    reciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
      (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ := by
  set g := psiGenStepResult K hq hπmax hπne0 f hf0 hf1 hf m (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m) with hg_def
  have hyψ : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
  have hyn : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (m + 1) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
  have hymem : (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure) ∈
      IntermediateField.adjoin K.carrier
        ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
          IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure)} : Set K.closure) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem
  haveI hyfd : FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier
      ({(↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
        IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure)} : Set K.closure)) := by
    rw [g.hcompat]; exact (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  have key := reciprocityMap_pred_eq_map_succ K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
    g.pt g.hψ g.hn g.hmem hyψ hyn hymem σ
  have hcongr := reciprocityMap_congr K hq hπmax hπne0 f hf0 hf1 hf (m + 1) (by omega)
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (m + 1 + 1) g.pt g.hn g.hmem π) :
      IntermediateField.adjoin K.carrier ({g.pt} : Set K.closure)) : K.closure)
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt g.hcompat hyψ hyn hymem
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hψ
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hn
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hmem σ
  exact key.trans hcongr

end ABC3.Found.PGC
