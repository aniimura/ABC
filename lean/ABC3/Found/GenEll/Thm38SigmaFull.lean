/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38SigmaZeta
import ABC3.Meta.Claim

/-!
# 第 1282 ブロック —— **`σ`・`ζ`・`π` の三つ組が揃う**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——第 1276 の入力が完成する

第 1276（`tateSigmaAct_unipotent`・`tateSigmaAct_exists_ne`）が要求するのは

| 仮説 | 意味 |
|---|---|
| `hζ : IsPrimitiveRoot ζ l` | `ζ` は原始 `l` 乗根 |
| `hπl : π^l = Q` | `π` は `Q` の `l` 乗根 |
| `hσζ : σU ζ = ζ` | `σ` は `ζ` を固定 |
| `hσπ : σU π = ζ π` | `σ` は `π` を `ζπ` に送る |

☆本ブロックは `l ∤ v(Q)` と「基礎体が `μ_l` を含む」から
**この 4 つを一度に出す**。

★★★`σ` が `ζ` を固定するのは、`ζ` が `l` 乗根なので基礎体に入るからである
——だから**大域の数体に `ζ_l` を添加しておく**（`ImageContainsSL2J` は
部分群の像で言えば十分なので、これは自由にできる）。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

/-- ★★★★★★★★★★★★★★★★★★★★
**`σ`・`ζ`・`π` の三つ組**——★**無条件**（第 1282）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`l ∤ v(Q)` から `σ`（第 1280）、`ζ`（第 1281）を取り、
基礎体が `μ_l` を含むことから `σ ζ = ζ` を得る。 -/
theorem exists_sigma_zeta_pi {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]
    [IsGalois K₀ M] (v₀ : K₀ˣ →* Multiplicative ℤ) {l : ℕ} (hl : l.Prime)
    (Q₀ : K₀ˣ) (hnd : ¬ ((l : ℤ) ∣ vAdd v₀ Q₀))
    (π : Mˣ) (hπ : π ^ l = Units.map (algebraMap K₀ M : K₀ →* M) Q₀)
    (hμ : ∀ z : Mˣ, z ^ l = 1 → ∃ y : K₀, algebraMap K₀ M y = ((z : M))) :
    ∃ (σ : M ≃ₐ[K₀] M) (ζ : Mˣ), IsPrimitiveRoot ((ζ : M)) l ∧
      Units.map σ.toAlgHom.toRingHom.toMonoidHom ζ = ζ ∧
      Units.map σ.toAlgHom.toRingHom.toMonoidHom π = ζ * π := by
  obtain ⟨σ, hσ⟩ := exists_algEquiv_move_pi v₀ Q₀ hnd π hπ
  obtain ⟨ζ, hζ, hσπ⟩ := sigma_pi_eq_zeta_mul hl Q₀ π hπ σ hσ
  refine ⟨σ, ζ, hζ, ?_, hσπ⟩
  have hζl : ζ ^ l = 1 := by
    refine Units.ext ?_
    push_cast
    simpa using hζ.pow_eq_one
  obtain ⟨y, hy⟩ := hμ ζ hζl
  refine Units.ext ?_
  show σ ((ζ : M)) = ((ζ : M))
  rw [← hy]
  exact σ.commutes y

/-! ## ★出典の紐付け(`.src`) -/

def exists_sigma_zeta_pi.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ・ζ・π の三つ組が揃う。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_sigma_zeta_pi.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_algEquiv_move_pi(第 1280、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_algEquiv_move_pi") 1,
    .citation "[ABC3]" "sigma_pi_eq_zeta_mul(第 1281、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sigma_pi_eq_zeta_mul") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1282）**——第 1276（`tateSigmaAct_unipotent`・" ++
       "`tateSigmaAct_exists_ne`）が要求する 4 つの仮説が一度に出る。" ++
       "☆`σ` が `ζ` を固定するのは `ζ` が `l` 乗根なので基礎体に入るからであり、" ++
       "そのために**大域の数体に `ζ_l` を添加しておく**" ++
       "（`ImageContainsSL2J` は部分群の像で言えば十分なので自由にできる）。") 2 ]

end ABC3.Found.GenEll
