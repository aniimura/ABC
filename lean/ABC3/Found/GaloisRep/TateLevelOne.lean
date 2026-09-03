/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BasisFree
import ABC3.Found.GaloisRep.TateWiring
import ABC3.Meta.Claim

/-!
# 第 1202 ブロック —— **`L̄` 上には位数 `l` の点があり、`T_l E` の第 1 層に持ち上がる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`HasLCyclicJ` を点へ運ぶ第 2 段

第 1201 で「`T_l E` の第 1 層の射影が `0` でなければ位数はちょうど `l`」を取った。
本ブロックはその**逆向き**と**存在**を取る:

| 定理 | 内容 | 材料 |
|---|---|---|
| `exists_point_addOrderOf_eq_prime` | ★★★`L̄` 上に位数 `l` の点がある | `torsion_card`（`#E[l] = l²`） |
| `exists_tateProj_one_eq` | ★★★`l`-捩れ点は `T_l E` の第 1 層に**持ち上がる** | `tateProj_surjective` ＋ `exists_smul_step` |

★★★これで `Lemma 3.5` を有限拡大の上で回す（第 1199）ための
**位数 `l` の点は必ず取れる**ことが言えた——`L̄` に取れば、
第 1195（`L(H)` は有限次）で有限拡大に落ちる。

☆残るのは `tateProj … 1` の**核**が `l • T_l E` であること
（`T_l E / l ≅ E[l]`）で、これは `HasLCyclicJ` の直線を点の直線へ
一意に対応させる段に要る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-! ## ★★★★★★★★★★★★位数 `l` の点は必ずある -/

/-- ★★★★★★★★★★★★**`L̄` 上には位数ちょうど `l` の点がある**——★**無条件**（第 1202）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`#E[l] = l² ≥ 4 > 1`（`torsion_card`）なので `0` でない `l`-捩れ点があり、
`l` が素数だからその位数はちょうど `l` である。

★★★これが `Lemma 3.5` を有限拡大の上で回す（第 1199）ための
**生成元の存在**である。 -/
theorem exists_point_addOrderOf_eq_prime [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    {l : ℕ} (hl : l.Prime) :
    ∃ P : (W.baseChange L).toAffine.Point, addOrderOf P = l := by
  have hΔ : (W.baseChange L).Δ ≠ 0 := (W.baseChange L).isUnit_Δ.ne_zero
  have hchar : ∀ k : ℕ, 1 ≤ k → ((k : L) ≠ 0) := fun k hk => Nat.cast_ne_zero.2 (by omega)
  have hcard : Nat.card {P : (W.baseChange L).toAffine.Point // l • P = 0} = l ^ 2 :=
    torsion_card (W.baseChange L) hΔ l hl.one_lt.le (fun k hk1 _ => hchar k hk1)
  have hne : ∃ P : (W.baseChange L).toAffine.Point, l • P = 0 ∧ P ≠ 0 := by
    by_contra hcon
    push_neg at hcon
    haveI : Subsingleton {P : (W.baseChange L).toAffine.Point // l • P = 0} :=
      ⟨fun a b => Subtype.ext (by rw [hcon a.1 a.2, hcon b.1 b.2])⟩
    haveI : Nonempty {P : (W.baseChange L).toAffine.Point // l • P = 0} :=
      ⟨⟨0, by simp⟩⟩
    have hone : Nat.card {P : (W.baseChange L).toAffine.Point // l • P = 0} = 1 :=
      Nat.card_eq_one_iff_unique.2 ⟨inferInstance, inferInstance⟩
    rw [hcard] at hone
    have h2 : 2 ≤ l := hl.two_le
    nlinarith [hone]
  obtain ⟨P, hP0, hPne⟩ := hne
  refine ⟨P, ?_⟩
  have hdvd := addOrderOf_dvd_of_nsmul_eq_zero hP0
  rcases hl.eq_one_or_self_of_dvd _ hdvd with h1 | h2
  · exact absurd (AddMonoid.addOrderOf_eq_one_iff.mp h1) hPne
  · exact h2

/-! ## ★★★★★★★★★★`l`-捩れ点は第 1 層に持ち上がる -/

/-- ★★★★★★★★★★**`l`-捩れ点は `T_l E` の第 1 層に持ち上がる**——★**無条件**（第 1202）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`tateProj_surjective`（第 205）の仮説 `hstep` は
`exists_smul_step`（第 209）で消える——標数 `0` の代数閉体では
`[l] : E[l^(m+1)] → E[l^m]` が全射だからである。

★★★これが第 1201 の `addOrderOf_tateProj_one` の**逆向き**であり、
`E[l]` の側の言明を `T_l E` の側へ運ぶ段である。 -/
theorem exists_tateProj_one_eq [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (l : ℕ) [Fact l.Prime]
    (P : (W.baseChange L).toAffine.Point) (hP : l • P = 0) :
    ∃ f : tateModule (W.baseChange L) l,
      ((tateProj (W.baseChange L) l 1 f : (W.baseChange L).toAffine.Point)) = P := by
  have hl : l.Prime := Fact.out
  have hΔ : (W.baseChange L).Δ ≠ 0 := (W.baseChange L).isUnit_Δ.ne_zero
  have hchar : ∀ k : ℕ, 1 ≤ k → ((k : L) ≠ 0) := fun k hk => Nat.cast_ne_zero.2 (by omega)
  have hstep : ∀ (m : ℕ) (T : (W.baseChange L).toAffine.Point), (l ^ m) • T = 0 →
      ∃ S : (W.baseChange L).toAffine.Point, (l ^ (m + 1)) • S = 0 ∧ l • S = T :=
    fun m T hT => exists_smul_step (W.baseChange L) hΔ l m hl.one_lt.le
      (fun k hk1 _ => hchar k hk1) hT
  have hP1 : (l ^ 1) • P = 0 := by simpa using hP
  obtain ⟨f, hf⟩ := tateProj_surjective (W.baseChange L) l 1 hstep ⟨P, hP1⟩
  exact ⟨f, congrArg Subtype.val hf⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_point_addOrderOf_eq_prime.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(L̄ 上には位数ちょうど l の点がある。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_point_addOrderOf_eq_prime.needs : List ProofObligation :=
  [ .citation "[ABC3]" "torsion_card(第 186、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.torsion_card") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1202）**——`Lemma 3.5` を有限拡大の上で回す" ++
       "（第 1199）ための**生成元の存在**である。" ++
       "☆`L̄` に取れば、第 1195（`L(H)` は有限次）で有限拡大に落ちる。") 2 ]

def exists_tateProj_one_eq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l-捩れ点は T_l E の第 1 層に持ち上がる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_tateProj_one_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateProj_surjective(第 205、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateProj_surjective") 1,
    .citation "[ABC3]" "exists_smul_step(第 209、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_smul_step") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1202）**——第 1201 の `addOrderOf_tateProj_one` の" ++
       "**逆向き**であり、`E[l]` の側の言明を `T_l E` の側へ運ぶ段である。" ++
       "☆残るのは `tateProj … 1` の**核**が `l • T_l E` であることで、" ++
       "これは `HasLCyclicJ` の直線を点の直線へ一意に対応させる段に要る。") 2 ]

end ABC3.Found.GaloisRep
