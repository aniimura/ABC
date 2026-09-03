/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalActVc
import ABC3.Meta.Claim

/-!
# 第 1289 ブロック —— **Tate モデルの上で Galois 作用は幂単かつ非自明**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——局所側の**まとめ**

第 1276（`tateSigmaAct` の幂単性・非自明性）と
第 1286（`tateSigmaAct` ＝ 実際の Galois 作用）を合わせるだけである。

☆受け取るのは:

| 仮説 | 出どころ |
|---|---|
| `TateSetup`・`Φ` | `mkTateSetup`・`dvrTatePhiAddEquiv`（在庫、無条件） |
| `ζ`（原始 `l` 乗根）・`π`（`π^l = Q`） | 第 1279（`E[l]` が載っていればよい） |
| `σ`（`σζ = ζ`・`σπ = ζπ`） | 第 1282（`l ∤ v(Q)` と `μ_l ⊆ K₀`） |

★★★これで**局所側は完成**である——残るのは
`E ⊗ K` と Tate モデルを結ぶ段（第 1288）と大域へ運ぶ段（第 1287・1270）だけ。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Tate モデルの上で Galois 作用は `E[l]` に幂単かつ非自明**——★**無条件**（第 1289）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1276 と第 1286 を合わせただけである。 -/
theorem galAct_unipotent_ne_of_tate (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+
      ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} [NeZero l] (hl1 : 1 < l) {ζ π : Kˣ}
    (hζ : IsPrimitiveRoot ((ζ : K)) l) (hπl : π ^ l = S.Q)
    (σR : R →+* R) (hσI : ∀ x ∈ I, σR x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σR r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (hUinj : Function.Injective σU)
    (hσζ : σU ζ = ζ) (hσπ : σU π = ζ * π)
    (hW : ((tateCurveAt S.q S.hq).map (algebraMap R K)).map σK
      = (tateCurveAt S.q S.hq).map (algebraMap R K)) :
    (∀ P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point, l • P = 0 →
        galAct σK _ hW (galAct σK _ hW P) + P = galAct σK _ hW P + galAct σK _ hW P) ∧
      (∃ P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point,
        l • P = 0 ∧ galAct σK _ hW P ≠ P) := by
  have heq := tateSigmaAct_eq_galAct S hΔ Φ hΦ σR hσI σK hcompat σU hσU hσq hσv hUinj hW
  constructor
  · intro P hP
    have h := tateSigmaAct_unipotent S hζ hπl Φ σU hσq hσζ hσπ P hP
    simp only [heq] at h
    exact h
  · obtain ⟨P, hP, hne⟩ := tateSigmaAct_exists_ne S hl1 hζ hπl Φ σU hσq hσπ
    refine ⟨P, hP, ?_⟩
    simp only [heq] at hne
    exact hne

/-! ## ★出典の紐付け(`.src`) -/

def galAct_unipotent_ne_of_tate.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate モデルの上で Galois 作用は E[l] に幂単かつ非自明。★無条件)",
    sectionId := "genell-thm-3-8" }

def galAct_unipotent_ne_of_tate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateSigmaAct_unipotent・tateSigmaAct_exists_ne(第 1276、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.tateSigmaAct_unipotent") 1,
    .citation "[ABC3]" "tateSigmaAct_eq_galAct(第 1286、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.tateSigmaAct_eq_galAct") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1289）**——**局所側の完成形**である。" ++
       "☆受け取るのは `TateSetup`・`Φ`（`mkTateSetup`・`dvrTatePhiAddEquiv`、在庫）、" ++
       "`ζ`・`π`（第 1279）、`σ`（第 1282）だけ。" ++
       "★残るのは `E ⊗ K` と Tate モデルを結ぶ段（第 1288）と" ++
       "大域へ運ぶ段（第 1287・1270）である。") 3 ]

end ABC3.Found.GenEll
