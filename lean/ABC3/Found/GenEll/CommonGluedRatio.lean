/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GluedGlobalRatio
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★試験元の族に対して**単一の指数**を取る —— 段 E3c の帳簿（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★これは何か

`§9-845` は**試験元 1 つ**につき指数 `n` と大域切断 `t` を作った。
★しかし `§9-843` の座標条件は**有限個の試験元すべて**について同時に要る
——座標族は 1 つの指数で揃っていなければならない。

★★本ファイルはその**指数揃え**である:

> 有限個の試験元 `g_k`（`k ∈ T`）に対し、**単一の `n`** と切断の族 `t_k` があって
> どの `k` でも `g_k = t_k / (s^{⊗(n+1)})`。

## ★★★機構 —— `exists_common_pow` を述語ごと当てる

★述語を「**分母も重なりも同時に満たす `a` がある**」と取るのが要点である:

    `P k n ≝ ∃ a, (∀ i, g_k·f_i^n = a_i) ∧ (∀ i j, a_i·u^n = a_j)`

★★これが `n` について単調であることは `res_pow_bump`（`§9-831`）と
`hagree_bump`（本ファイル、`§9-840` の本体）から出る
——`a_i ↦ a_i·f_i^m` で両方が同時に持ち上がる。
★★★あとは `exists_common_pow`（`§9-826`、`Finset.sup`）に渡すだけである。

## ★これで段 E3c の帳簿が閉じる

★★残るのは「`t_k` たちと `s^{⊗(n+1)}` を `Fin (N+1)` に並べる」という
**添字の付け替えだけ**であり、そこに数学は無い。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★重なりの一致条件は底上げで保たれる -/

/-- ★★★**重なりの一致条件は指数の底上げで保たれる**（`§9-840` の本体を補題にしたもの）。

★`f_j| = f_i|·u` なので `u^m` が両辺に現れてちょうど相殺する。 -/
theorem hagree_bump {ι : Type} (M : X.PresheafOfModules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (n m : ℕ) (a : ∀ i, (Γ(X, U i) : Type))
    (hagree : ∀ i j,
      X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) :
    ∀ i j,
      X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op
          (a i * (trivValue M (U i) (e i) s) ^ m)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ (n + m)
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op
          (a j * (trivValue M (U j) (e j) s) ^ m) := by
  intro i j
  simp only [map_mul, map_pow]
  rw [trivValue_res_transUnit M (U i) (U j) (e i) (e j) s]
  refine overlap_bump _ _ _ _ n m ?_
  rw [hagree i j]

/-! ## ★★★★★★★★試験元の族に対する単一の指数 -/

open scoped Classical in
/-- ★★★★★★★★**試験元の族に対して単一の指数を取る**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

有限個の試験元 `g_k`（`k ∈ T`）に対し、★**単一の `n`** と分子の族 `a_k` があって
どの `k` でも「分母が払えている」かつ「重なりで一致する」。

★★機構は `exists_common_pow`（`§9-826`）に**述語ごと**当てることである
——述語を「分母も重なりも同時に満たす `a` がある」と取れば、
単調性は `res_pow_bump` と `hagree_bump` から出る。 -/
theorem exists_common_exponent_family {ι κ : Type} [Fintype ι] (T : Finset κ)
    (M : X.PresheafOfModules) (U : ι → X.Opens)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (g : κ → (Γ(X, nonVanishing M s) : Type)) :
    ∃ (n : ℕ) (a : κ → ∀ i, (Γ(X, U i) : Type)), ∀ k ∈ T,
      (∀ i, X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op (g k)
          * (algebraMap (Γ(X, U i) : Type)
              (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type)
              (trivValue M (U i) (e i) s)) ^ n
        = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type) (a k i)) ∧
      (∀ i j, X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a k i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a k j)) := by
  have key : ∃ n : ℕ, ∀ k ∈ T, ∃ a : ∀ i, (Γ(X, U i) : Type),
      (∀ i, X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op (g k)
          * (algebraMap (Γ(X, U i) : Type)
              (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type)
              (trivValue M (U i) (e i) s)) ^ n
        = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type) (a i)) ∧
      (∀ i j, X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) := by
    refine exists_common_pow T _ ?_ ?_
    · rintro k n₁ n₂ hle ⟨a, hden, hag⟩
      obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le hle
      exact ⟨fun i => a i * (trivValue M (U i) (e i) s) ^ m,
        fun i => res_pow_bump (U i) _ _ (a i) n₁ m (hden i),
        hagree_bump M U e s n₁ m a hag⟩
    · intro k _
      obtain ⟨n, a, hden, hag⟩ := exists_common_exponent Finset.univ M U
        (fun i _ => hU i) (fun i _ j _ => hUij i j) e s
        (fun i => X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op (g k))
        (fun i _ j _ => by rw [res_trans, res_trans])
      exact ⟨n, a, fun i => hden i (Finset.mem_univ i),
        fun i j => hag i (Finset.mem_univ i) j (Finset.mem_univ j)⟩
  obtain ⟨n, hn⟩ := key
  choose! a ha using hn
  exact ⟨n, a, fun k hk => ha k hk⟩

/-! ## ★★★★★★★★★分子と分母が揃っていれば大域の比である -/

set_option maxHeartbeats 4000000 in
/-- ★★★★★★★★★**分子と分母が揃っていれば大域の比である**（`§9-845` の尾部を補題に）。

★`§9-845` はこの補題の系である。族に対しても同じ尾部を使うので、名前を付けて括り出した。 -/
theorem globalRatio_of_den_agree {ι : Type} (M : X.PresheafOfModules)
    (hM : IsLocallyTrivial X M) (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (n : ℕ) (g : (Γ(X, nonVanishing M s) : Type))
    (a : ∀ i, (Γ(X, U i) : Type))
    (t : (((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val.obj
      (op (⊤ : X.Opens)) : Type))
    (hden : ∀ i, X.presheaf.map (homOfLE (basicOpen_trivValue_le M (U i) (e i) s)).op g
          * (algebraMap (Γ(X, U i) : Type)
              (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type)
              (trivValue M (U i) (e i) s)) ^ n
        = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) s)) : Type) (a i))
    (ht : ∀ i, trivValue ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val (U i)
        (sheafifyTriv (presheafTensorPow M (n + 1)) (tensorPowTriv (e i) (n + 1))) t
      = a i * (trivValue M (U i) (e i) s) ^ 1) :
    X.presheaf.map (homOfLE (le_of_eq (nonVanishing_unit_secPow M hM s n))).op g
      = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
          (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
          t (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
            (secPow M s (n + 1))) := by
  refine globalRatio_unique_of_cover _ _ t _ U
    (fun i => sheafifyTriv (presheafTensorPow M (n + 1)) (tensorPowTriv (e i) (n + 1)))
    (by rw [← inf_iSup_eq, hcov, inf_top_eq]) _ (fun i => ?_)
  refine sectionRatio_unique _ (U i) _ t _ _ ?_
  have hOpen := nonVanishing_inf_unit_secPow M (U i) (e i) s n
  have hd := congrArg (fun z => X.presheaf.map (homOfLE (le_of_eq hOpen)).op z) (hden i)
  simp only [map_mul, map_pow, algebraMap_basicOpen_eq_res] at hd
  rw [res_trans, res_trans, res_trans] at hd
  rw [trivValue_unit_secPow, ht i, pow_one, map_mul, map_pow, res_trans,
    pow_succ, ← mul_assoc, hd]

/-! ## ★★★★★★★★★★試験元の族が同時に大域の比になる -/

open scoped Classical in
/-- ★★★★★★★★★★**有限個の試験元が同時に大域の比になる** —— 段 E3c の帳簿。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

有限個の試験元 `g_k`（`k ∈ T`）に対し、★**単一の `n`** と切断の族 `t_k` があって

    `g_k = t_k / (s^{⊗(n+1)})`   （大域の比として、どの `k` でも同じ `n`）

★★これで `§9-843` の座標条件が**族ごと**満たせる
——座標族に `s^{⊗(n+1)}` と `t_k` たちを並べればよい。
★★★残るのは `Fin (N+1)` への**添字の付け替えだけ**であり、そこに数学は無い。 -/
theorem exists_common_glued_globalRatio {ι κ : Type} [Fintype ι] (T : Finset κ)
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (s : (M.obj (op ⊤) : Type)) (g : κ → (Γ(X, nonVanishing M s) : Type)) :
    ∃ (n : ℕ) (t : κ → (((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val.obj
        (op (⊤ : X.Opens)) : Type)), ∀ k ∈ T,
      X.presheaf.map (homOfLE (le_of_eq (nonVanishing_unit_secPow M hM s n))).op (g k)
        = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
            (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
            (t k) (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
              (secPow M s (n + 1))) := by
  obtain ⟨n, a, hka⟩ := exists_common_exponent_family T M U hU hUij e s g
  have hT : ∀ k ∈ T, ∃ t : (((sheafifyFunctor X).obj
      (presheafTensorPow M (n + 1))).val.obj (op (⊤ : X.Opens)) : Type),
      ∀ i, trivValue ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val (U i)
          (sheafifyTriv (presheafTensorPow M (n + 1)) (tensorPowTriv (e i) (n + 1))) t
        = a k i * (trivValue M (U i) (e i) s) ^ 1 :=
    fun k hk => exists_glue_bump_trivValue M U hcov e s n 1 (a k) (hka k hk).2
  choose! t ht using hT
  exact ⟨n, t, fun k hk =>
    globalRatio_of_den_agree M hM U hcov e s n (g k) (a k) (t k) (hka k hk).1 (ht k hk)⟩

/-! ## ★出典の紐付け(`.src`) -/

def hagree_bump.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの一致条件は指数の底上げで保たれる)",
    sectionId := "genell-prop-1-4" }

def exists_common_exponent_family.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(試験元の族に対して単一の指数を取る)",
    sectionId := "genell-prop-1-4" }

def globalRatio_of_den_agree.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分子と分母が揃っていれば大域の比である)",
    sectionId := "genell-prop-1-4" }

def exists_common_glued_globalRatio.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(有限個の試験元が同時に大域の比になる)",
    sectionId := "genell-prop-1-4" }

def exists_common_glued_globalRatio.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_common_pow(指数を最大値で揃える、§9-826)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_pow") 2,
    .citation "[ABC3]" "exists_common_exponent(§9-831) / exists_glue_bump_trivValue(§9-840)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_exponent") 2,
    .citation "[ABC3]" "globalRatio_unique_of_cover(§9-844)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_unique_of_cover") 2,
    .implicitStep
      ("★述語を「**分母も重なりも同時に満たす a がある**」と取るのが要点である" ++
       "——単調性は res_pow_bump(§9-831)と hagree_bump(本ファイル)から出る" ++
       "(a_i ↦ a_i·f_i^m で両方が同時に持ち上がる)") 4,
    .implicitStep
      ("★★★残るのは Fin (N+1) への**添字の付け替えだけ**であり、そこに数学は無い" ++
       "——座標族に s^{⊗(n+1)} と t_k たちを並べる帳簿である") 4 ]

end ABC3.Found.GenEll
