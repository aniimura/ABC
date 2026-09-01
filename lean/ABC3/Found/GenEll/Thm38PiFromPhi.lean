/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38TorsionNotMu
import ABC3.Meta.Claim

/-!
# 第 1279 ブロック —— **`E[l]` が載っていれば `π` は作れる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——第 1277・1278 の合流

☆Tate 一意化 `Φ : Additive (Kˣ ⧸ ⟨Q⟩) ≃+ P` があり、
`P` の `l`-捩れが `l` 個より多ければ、**`Q` は `K` の中で `l` 乗になる**。

★機構は 2 段:

| 段 | 内容 | 第 |
|---|---|---|
| 1 | `l`-捩れが `l` 個より多ければ `μ_l` の像からはみ出す類がある | 1278 |
| 2 | はみ出す類の代表 `x` は `x^l = Q^k`（`l ∤ k`）を満たし、Bezout で `π` が出る | 1277 |

★★★これで **`π` を作るために体を拡げる必要は無い**——
`E[l]` が `K` の上に載ってさえいればよい。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] {I : Ideal R} {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★★★★★★★★★★★★★
**`E[l]` が `l` 個より多く載っていれば `Q` は `l` 乗**——★**無条件**（第 1279）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆消費側では `T` を `E[l]`（個数 `l²`、`torsion_card`）から取る。 -/
theorem exists_pi_of_phi (S : TateSetup R I K) {l : ℕ} (hl : l.Prime)
    {ζ : Kˣ} (hζ : IsPrimitiveRoot ((ζ : K)) l)
    {P : Type*} [AddCommGroup P] (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+ P)
    (T : Finset P) (hT : ∀ p ∈ T, l • p = 0) (hcard : l < T.card) :
    ∃ π : Kˣ, π ^ l = S.Q := by
  classical
  haveI : NeZero l := ⟨hl.ne_zero⟩
  have hinj : Function.Injective
      (fun p : P => Additive.toMul (Φ.symm p)) := by
    intro a b hab
    exact Φ.symm.injective (Additive.toMul.injective hab)
  have hcard' : (T.image (fun p : P => Additive.toMul (Φ.symm p))).card = T.card :=
    Finset.card_image_of_injective T hinj
  have hT'l : ∀ c ∈ T.image (fun p : P => Additive.toMul (Φ.symm p)), c ^ l = 1 := by
    intro c hc
    obtain ⟨p, hp, rfl⟩ := Finset.mem_image.1 hc
    have h0 : l • (Φ.symm p) = 0 := by
      rw [← map_nsmul, hT p hp, map_zero]
    have h1 := congrArg Additive.toMul h0
    simpa using h1
  obtain ⟨c, hc1, hc2⟩ := exists_torsion_not_mu hl.pos S.Q ζ
    (units_pow_eq_one_of_isPrimitiveRoot hζ)
    (T.image (fun p : P => Additive.toMul (Φ.symm p))) hT'l (by omega)
  exact exists_pow_eq_of_torsion_not_mu hl S.Q ζ
    (fun y hy => exists_zpow_of_pow_eq_one hζ y hy) c hc1 hc2

/-! ## ★出典の紐付け(`.src`) -/

def exists_pi_of_phi.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(E[l] が l 個より多く載っていれば Q は l 乗。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_pi_of_phi.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_torsion_not_mu(第 1278、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_torsion_not_mu") 1,
    .citation "[ABC3]" "exists_pow_eq_of_torsion_not_mu(第 1277、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pow_eq_of_torsion_not_mu") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1279）**——**`π` を作るために体を拡げる必要は無い**。" ++
       "☆`E[l]` が `K` の上に載ってさえいればよく、それは大域の数体を" ++
       "`L′ ≔ L(ζ_l, E[l], √d)` に取り替えれば満たされる。" ++
       "★★★これが第 1274 の B1（完備 DVR の Eisenstein 拡大）を道から消した理由である。") 2 ]

end ABC3.Found.GenEll
