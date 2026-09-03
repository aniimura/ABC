/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.TateGalActUnip
import ABC3.Found.GaloisRep.TateVarChange
import ABC3.Meta.Claim

/-!
# 第 1290 ブロック —— **2 条件は加法同型で移る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——Tate モデルから `E` へ戻す段

第 1289 で「Tate モデルの上で Galois 作用は幂単かつ非自明」を取った。
★Tate モデルは `E ⊗ K` の変数変換なので、
**加法同型で 2 条件を戻せば `E ⊗ K` の話になる**。

☆本ブロックはその**抽象な形**（`e : A ≃+ B` と作用の絡み合い）と、
変数変換への当てはめである。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★
**幂単性・非自明性は加法同型で移る**——★**無条件**（第 1290）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`e (τA a) = τB (e a)` なら、`B` 側の 2 条件が `A` 側に降りる。 -/
theorem unipotent_ne_transfer {A B : Type*} [AddCommGroup A] [AddCommGroup B]
    (e : A ≃+ B) (τA : A →+ A) (τB : B →+ B) (hint : ∀ a, e (τA a) = τB (e a)) (n : ℕ)
    (h1 : ∀ b : B, n • b = 0 → τB (τB b) + b = τB b + τB b)
    (h2 : ∃ b : B, n • b = 0 ∧ τB b ≠ b) :
    (∀ a : A, n • a = 0 → τA (τA a) + a = τA a + τA a) ∧
      (∃ a : A, n • a = 0 ∧ τA a ≠ a) := by
  constructor
  · intro a ha
    refine e.injective ?_
    have hA : e (τA (τA a)) = τB (τB (e a)) := by rw [hint, hint]
    have hB : e (τA a) = τB (e a) := hint a
    rw [map_add, map_add, hA, hB]
    exact h1 (e a) (by rw [← map_nsmul, ha, map_zero])
  · obtain ⟨b, hb, hne⟩ := h2
    refine ⟨e.symm b, ?_, ?_⟩
    · rw [← map_nsmul, hb, map_zero]
    · intro hcon
      refine hne ?_
      have h3 := congrArg e hcon
      rw [hint, e.apply_symm_apply] at h3
      exact h3

/-- ★★★★★★**変数変換の加法同型は `GenEll.vcPoint` そのもの**——★**無条件**（第 1290）。

☆`GaloisRep.vcPoint W C` と `GenEll.vcPoint C W` は引数の順が違うだけである。 -/
theorem vcPointAddEquiv_apply {F : Type} [Field F] (W : WeierstrassCurve F)
    (C : VariableChange F) (P : W.toAffine.Point) :
    ABC3.Found.GaloisRep.vcPointAddEquiv W C P = vcPoint C W P := rfl

/-- ★★★★★★★★★★★★★★★★
**Tate モデルの 2 条件は `E ⊗ K` に戻る**——★**無条件**（第 1290）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆絡み合いは第 1288（`galAct` は変数変換と可換）である。 -/
theorem unipotent_ne_of_variableChange {K : Type} [Field K] (σK : K →+* K)
    (W : WeierstrassCurve K) (C : VariableChange K)
    (hC : C.map σK = C) (hW : W.map σK = W) (n : ℕ)
    (h1 : ∀ P : (C • W).toAffine.Point, n • P = 0 →
      galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW])
          (galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW]) P) + P
        = galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW]) P
          + galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW]) P)
    (h2 : ∃ P : (C • W).toAffine.Point, n • P = 0 ∧
      galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW]) P ≠ P) :
    (∀ P : W.toAffine.Point, n • P = 0 →
        galAct σK W hW (galAct σK W hW P) + P = galAct σK W hW P + galAct σK W hW P) ∧
      (∃ P : W.toAffine.Point, n • P = 0 ∧ galAct σK W hW P ≠ P) :=
  unipotent_ne_transfer (ABC3.Found.GaloisRep.vcPointAddEquiv W C) (galAct σK W hW)
    (galAct σK (C • W) (by rw [← WeierstrassCurve.map_variableChange, hC, hW]))
    (fun P => by
      rw [vcPointAddEquiv_apply, vcPointAddEquiv_apply]
      exact (galAct_vcPoint σK W C hC hW P).symm)
    n h1 h2

/-! ## ★出典の紐付け(`.src`) -/

def unipotent_ne_transfer.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(幂単性・非自明性は加法同型で移る。★無条件)",
    sectionId := "genell-thm-3-8" }

def vcPointAddEquiv_apply.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(変数変換の加法同型は GenEll.vcPoint そのもの。★無条件)",
    sectionId := "genell-thm-3-8" }

def unipotent_ne_of_variableChange.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate モデルの 2 条件は E ⊗ K に戻る。★無条件)",
    sectionId := "genell-thm-3-8" }

def unipotent_ne_of_variableChange.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galAct_vcPoint(第 1288、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galAct_vcPoint") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1290）**——Tate モデルから `E` へ戻す段である。" ++
       "☆絡み合いは第 1288（`galAct` は変数変換と可換）。" ++
       "★これで局所側は `E ⊗ K` の言葉になった。") 2 ]

end ABC3.Found.GenEll
