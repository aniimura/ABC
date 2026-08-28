/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.DegFinBaseChange
import ABC3.Found.Arakelov.ArchDegBaseChange
import ABC3.Found.Arakelov.DegAPicM
import ABC3.Found.Arakelov.HeightUncond
import ABC3.Found.Arakelov.PullbackFunctorial
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★**底変換が閉じた** `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` と `ht_M̄` の整合（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L)

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

原文が

> `deg_K(L|_{Spec(𝓞_K)}) = deg_F(L)` ... In particular, ... it makes sense to define
> `ht_M̄(x) ≝ deg_F(x_F^* M̄)` — where `x_F` is **any** morphism that gives rise to `x`.

と書く段である。★**`ht_M̄` が `x_F` の取り方に依らない**ことがここで出る。

## ★★★★★機構（4 段）

    degFinPre_baseChange : 有限側 `degFin_K = [K:F]·degFin_F`（`§9-798`）
    archDeg_baseChange   : アルキメデス側 `archDeg_K = archDeg_F`（`§9-796`）
    ★正規化 `1/[K:ℚ] = 1/([F:ℚ]·[K:F])` が `[K:F]` を打ち消す
    APicMPullback_comp   : 引き戻しの関手性（`§9-747`）

★★`pullSec s ≠ 0` は**位数の計算から出る**（`pullSecMod_ne_zero`）——
もし `0` なら商が全体になり位数 `0`（無限）になるが、
`§9-798` はそれが正の数の `[K:F]` 乗だと言うからである。
★★★可逆加群は無限（`infinite_of_invertible`）——係数環が無限だからである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace NumberField
open ABC3.Found.GenEll

/-- ★**可逆加群は無限である**（係数環が無限な整域なら）。 -/
theorem infinite_of_invertible (S : Type) [CommRing S] [IsDomain S] [Infinite S]
    (A : Type) [AddCommGroup A] [Module S A] [Module.Invertible S A] : Infinite A := by
  haveI : Nontrivial A := nontrivial_of_invertible S A
  obtain ⟨m, hm⟩ := exists_ne (0 : A)
  exact Infinite.of_injective (fun r => (LinearMap.toSpanSingleton S A m) r)
    (toSpanSingleton_injective_invertible S A m hm)

set_option maxHeartbeats 2000000 in
/-- ★★★★★**非零切断の引き戻しは非零である**。

★機構は**位数の計算**である——もし `0` なら商が全体になり位数 `0`（無限）だが、
`§9-798` の `card_quotient_baseChange` はそれが正の数の `[K:F]` 乗だと言う。 -/
theorem pullSecMod_ne_zero (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)) (hs : s ≠ 0) :
    pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s ≠ 0 := by
  intro h
  have hc := card_quotient_baseChange F K L s hs
  rw [h, Submodule.span_zero_singleton] at hc
  haveI hinvF := invertible_gammaModPre (CommRingCat.of (𝓞 F)) L
  haveI hQ : Finite ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
      ⧸ (((CommRingCat.of (𝓞 F)) : Type) ∙ s)) :=
    finite_quotient_span_invertible ((CommRingCat.of (𝓞 F)) : Type)
      (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
      (fun r hr => finite_quotient_span F r hr) s hs
  have hpos0 : 0 < Nat.card ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
      ⧸ (((CommRingCat.of (𝓞 F)) : Type) ∙ s)) := Nat.card_pos
  have hpos : 0 < (Nat.card ((gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)
      ⧸ ((Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type) ∙ s)))
      ^ (Module.finrank F K) := by
    rw [card_quotient_span_gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf s]
    exact pow_pos hpos0 _
  rw [← hc] at hpos
  have e0 := Submodule.quotEquivOfEqBot
    (⊥ : Submodule (Γ(Spec (CommRingCat.of (𝓞 K)),
      (⊤ : (Spec (CommRingCat.of (𝓞 K))).Opens)) : Type)
      ((gammaModPre (CommRingCat.of (𝓞 K))
        ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))))).obj
          L.carrier.sheaf) : Type))) rfl
  rw [Nat.card_congr e0.toEquiv] at hpos
  haveI hinvK : Module.Invertible ((CommRingCat.of (𝓞 K)) : Type)
      ((gammaModPre (CommRingCat.of (𝓞 K))
        ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))))).obj
          L.carrier.sheaf) : Type)) :=
    invertible_gammaModPre (CommRingCat.of (𝓞 K))
      (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
  haveI : Infinite ((gammaModPre (CommRingCat.of (𝓞 K))
      ((pullbackPre (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))))).obj
        L.carrier.sheaf) : Type)) :=
    infinite_of_invertible ((CommRingCat.of (𝓞 K)) : Type) _
  rw [Nat.card_eq_zero_of_infinite] at hpos
  exact absurd hpos (lt_irrefl 0)

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★**`deg_F` は底変換で不変である**（切断を指定した形）。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L) -/
theorem degArithPre_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (gammaModPre (CommRingCat.of (𝓞 F)) L.carrier.sheaf : Type)) (hs : s ≠ 0) :
    degArithPre K (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
        (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s)
      = degArithPre F L s := by
  have hfin := degFinPre_baseChange F K L s hs
  have harch : archDeg K
      (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L).carrier
      (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s) = archDeg F L.carrier s :=
    archDeg_baseChange F K L.carrier s
  have hFK : (Module.finrank ℚ F) * (Module.finrank F K) = Module.finrank ℚ K :=
    Module.finrank_mul_finrank ℚ F K
  have hFpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  have hKFpos : (0 : ℝ) < (Module.finrank F K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := F) (M := K)
  show degFinPre (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
      (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s) / (Module.finrank ℚ K : ℝ)
    + archDeg K (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L).carrier
      (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf s) = _
  rw [hfin, harch]
  show _ = degFinPre L s / (Module.finrank ℚ F : ℝ) + archDeg F L.carrier s
  rw [← hFK]
  push_cast
  field_simp

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**`deg_K(L|_{Spec 𝓞_K}) = deg_F(L)`**——原文そのもの。

原文 (GenEll p.4):
> degK(L|Spec(OK)) = degF (L) -/
theorem degAInv_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    degAInv K (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
      = degAInv F L := by
  have hs := chosenSec_ne_zero F L
  have hne := pullSecMod_ne_zero F K L (chosenSec F L) hs
  rw [degAInv_eq K (AInv.pullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) L)
      (pullSecMod (𝓞 F) (𝓞 K) L.carrier.sheaf (chosenSec F L)) hne,
    degAInv_eq F L (chosenSec F L) hs]
  exact degArithPre_baseChange F K L (chosenSec F L) hs

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★**類の水準での底変換不変性**。 -/
theorem degAPicM_baseChange (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (a : APicM (Spec (CommRingCat.of (𝓞 F)))) :
    degAPicM K (APicMPullback (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K)))) a)
      = degAPicM F a := by
  induction a using Quotient.ind with
  | _ L => exact degAInv_baseChange F K L

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**高さは `x_F` の取り方に依らない**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★これが原文の「**any** morphism that gives rise to `x`」の中身である
——`x_F` を体の拡大で取り替えても `ht_M̄` は変わらない。
★★★機構は `§9-747` の関手性（`APicMPullback_comp`）と `degAPicM_baseChange` の合成だけである。 -/
theorem htMetricU_baseChange {X : Scheme.{0}} (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (M : AInv X)
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU K M (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF)
      = htMetricU F M xF := by
  show degAPicM K (APicMPullback
    (Spec.map (CommRingCat.ofHom (algebraMap (𝓞 F) (𝓞 K))) ≫ xF) (APicM.mk M)) = _
  rw [← APicMPullback_comp]
  exact degAPicM_baseChange F K _

/-! ### ★出典の紐付け(`.src`) -/

def degAInv_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_K(L|Spec(𝓞_K)) = deg_F(L)——正規化次数の底変換不変性)",
    sectionId := "genell-def-1-1-ii" }

def htMetricU_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(高さは x_F の取り方に依らない)",
    sectionId := "genell-def-1-1-ii" }

def htMetricU_baseChange.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degFinPre_baseChange(有限側、§9-798)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degFinPre_baseChange") 4,
    .citation "[ABC3]" "archDeg_baseChange(アルキメデス側、§9-796)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.archDeg_baseChange") 4,
    .citation "[ABC3]" "APicMPullback_comp(引き戻しの関手性、§9-747)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicMPullback_comp") 3,
    .implicitStep
      ("★原文は deg_F を ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)([Szp] Prop 1.1)で移して" ++
       "定義するが、本構成はその同型を**経由しない**(逸脱の記録、§9-782)。" ++
       "★★消費側(Definition 1.2 の ht)が要るのは deg_F が準同型で底変換で不変なことだけなので、" ++
       "後続の証明に影響は出ない") 4 ]

end ABC3.Found.Arakelov
