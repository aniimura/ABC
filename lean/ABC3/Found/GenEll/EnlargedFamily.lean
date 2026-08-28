/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CommonGluedRatioMulti
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★族を拡大する —— 段 E3 の組み立て（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★★これは何か —— 「分子を足す」を 1 本の補題に

段 E3（チャート写像 `A⁰_{x_i} → Γ(X, X_{s_i})` の全射性）の道具は

* `exists_finset_surjective_globalAwayHom`（`§9-847`）——
  「有限個の試験元 `T_i` が全部 `s_j/s_i` の形なら全射」
* `exists_common_glued_globalRatio_multi`（`§9-917`）——
  「単一の `n` と分子 `t_k` で `g_k = t_k / s_{d k}^{⊗(n+1)}`」

★★★★本ファイルはこの 2 つを繋ぐ**族の拡大**を 1 本にする:

    分母 `{s_i}` と試験元 `{T_i}` ⟹ 拡大した族 `s' : Fin (N'+1) → Γ(M^{⊗(n+1)})`
      * `s'(ρ i) = (s_i)^{⊗(n+1)}`（分母は残る）
      * どの `g ∈ T_i` も `s'_j / s'(ρ i)` の形になる

## ★★★機構 —— 直和して `Fin` に並べ直す

★添字を `Fin (N+1) ⊕ (Σ i, T_i)` に取り、`exists_fin_reindex` で `Fin (N'+1)` に移す。
★★`Fin (N+1)` が空でないので `Fin (N+1) ⊕ κ` も空でなく、`N'+1 = card` が取れる。

## ★★これで何が保たれるか

★★★**非消失軌跡は変わらない**——`nonVanishing_unit_secPow` により
`X_{(s_i)^{⊗(n+1)}} = X_{s_i}` である。したがって

* **被覆**（`⨆_i X_{s_i} = ⊤`）はそのまま
* **アフィン性**（`§9-915`）もそのまま
* ★分子 `t_k` の非消失軌跡はアフィンとは限らないが、
  `§9-916`（埋め込み性は部分族で足りる）があるので構わない

## ★配管の記録

★★★★`globalRatio _ _ (s' j) (s' (ρ i))` と書くと**型が合わない**
——`globalRatio` の型は分母 `s' (ρ i)` に依存するのに、
左辺の型は `nonVanishing_unit_secPow` を通して `(s_i)^{⊗(n+1)}` の形で書かれているからである。
★分母は `s' (ρ i)` ではなく `(s_i)^{⊗(n+1)}` と**明示的に**書く（等しいことは第 1 結論）。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★有限型の族を `Fin (N+1)` に並べ直す -/

/-- ★**空でない有限型の族は `Fin (N+1)` に並べ直せる**。 -/
theorem exists_fin_reindex {ι A : Type} [Fintype ι] [Nonempty ι] (s : ι → A) :
    ∃ (N : ℕ) (s' : Fin (N + 1) → A) (ρ : ι → Fin (N + 1)),
      (∀ i, s' (ρ i) = s i) ∧ (∀ j, ∃ i, s' j = s i) := by
  classical
  obtain ⟨N, hN⟩ : ∃ N, Fintype.card ι = N + 1 :=
    ⟨Fintype.card ι - 1, by have := Fintype.card_pos (α := ι); omega⟩
  let e : ι ≃ Fin (N + 1) := (Fintype.equivFin ι).trans (finCongr hN)
  refine ⟨N, fun j => s (e.symm j), e, fun i => by simp [e], fun j => ⟨e.symm j, rfl⟩⟩

/-! ## ★★★★★★★★★★★★★族の拡大 -/

set_option maxHeartbeats 2000000 in
open scoped Classical in
/-- ★★★★★★★★★★★★★**族を拡大して試験元を比にする** —— 段 E3 の組み立て。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

分母の族 `{s_i}`（`Fin (N+1)` で添字づけ）と各 `i` の試験元 `T_i` に対し、
★**単一の `n`** と拡大した族 `s' : Fin (N'+1) → Γ(X, sheafify(M^{⊗(n+1)}))` があって

* `s' (ρ i) = (s_i)^{⊗(n+1)}` ——**分母はそのまま残る**
* どの `g ∈ T_i` も `s'_j / (s_i)^{⊗(n+1)}` の形になる

★★これを `§9-847` の `exists_finset_surjective_globalAwayHom` に渡すと
**チャート写像の全射性**が出る。 -/
theorem exists_enlarged_family {ι : Type} [Fintype ι]
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    {N : ℕ} (s : Fin (N + 1) → (M.obj (op ⊤) : Type))
    (Tfam : ∀ i : Fin (N + 1), Finset ((Γ(X, nonVanishing M (s i)) : Type))) :
    ∃ (n N' : ℕ)
      (s' : Fin (N' + 1) → (((sheafifyFunctor X).obj
        (presheafTensorPow M (n + 1))).val.obj (op (⊤ : X.Opens)) : Type))
      (ρ : Fin (N + 1) → Fin (N' + 1)),
      (∀ i, s' (ρ i) = ((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
          (secPow M (s i) (n + 1))) ∧
      (∀ i, ∀ g ∈ Tfam i, ∃ j : Fin (N' + 1),
        X.presheaf.map (homOfLE (le_of_eq (nonVanishing_unit_secPow M hM (s i) n))).op g
          = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
              (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
              (s' j) (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
                (secPow M (s i) (n + 1)))) := by
  set κ : Type := Σ i : Fin (N + 1), {g // g ∈ Tfam i} with hκ
  obtain ⟨n, t, ht⟩ := exists_common_glued_globalRatio_multi (Finset.univ : Finset κ)
    M hM U hcov hU hUij e s (fun k => k.1) (fun k => k.2.1)
  obtain ⟨N', s', ρ', hρ'1, -⟩ := exists_fin_reindex
    (Sum.elim (fun i : Fin (N + 1) =>
      ((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
        (secPow M (s i) (n + 1))) t)
  refine ⟨n, N', s', fun i => ρ' (Sum.inl i), fun i => hρ'1 (Sum.inl i), ?_⟩
  intro i g hg
  refine ⟨ρ' (Sum.inr ⟨i, ⟨g, hg⟩⟩), ?_⟩
  rw [hρ'1 (Sum.inr ⟨i, ⟨g, hg⟩⟩)]
  exact ht ⟨i, ⟨g, hg⟩⟩ (Finset.mem_univ _)

/-! ## ★出典の紐付け(`.src`) -/

def exists_fin_reindex.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(空でない有限型の族は Fin (N+1) に並べ直せる)",
    sectionId := "genell-prop-1-4" }

def exists_enlarged_family.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(族を拡大して試験元を比にする)",
    sectionId := "genell-prop-1-4" }

def exists_enlarged_family.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_common_glued_globalRatio_multi(分母が複数でも単一の指数、§9-917)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_glued_globalRatio_multi") 3,
    .citation "[ABC3]" "nonVanishing_unit_secPow(冪を上げても非消失軌跡は変わらない)"
      (.inProject "ABC3" "ABC3.Found.GenEll.nonVanishing_unit_secPow") 1,
    .implicitStep
      ("★★★非消失軌跡が変わらないので、被覆(⨆_i X_{s_i} = ⊤)も" ++
       "アフィン性(§9-915)もそのまま保たれる。" ++
       "★分子 t_k の非消失軌跡はアフィンとは限らないが、" ++
       "§9-916(埋め込み性は部分族で足りる)があるので構わない") 3,
    .implicitStep
      ("★★★★配管: globalRatio _ _ (s' j) (s' (ρ i)) と書くと**型が合わない**" ++
       "——globalRatio の型は分母 s' (ρ i) に依存するのに、" ++
       "左辺の型は nonVanishing_unit_secPow を通して (s_i)^{⊗(n+1)} の形で" ++
       "書かれているからである。分母は明示的に書く(等しいことは第 1 結論)") 2,
    .implicitStep
      ("★★次は これを §9-847 の exists_finset_surjective_globalAwayHom に渡して" ++
       "チャート写像の全射性を出し、§9-916 に渡す段である") 4 ]

end ABC3.Found.GenEll
