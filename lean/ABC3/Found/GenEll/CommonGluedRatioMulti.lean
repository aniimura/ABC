/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.CommonGluedRatio
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★分母が**複数**でも単一の指数が取れる —— 段 E3 の要（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★★★★★これは何か —— 段 E3 の組み立てを可能にする

`§9-CommonGluedRatio` の `exists_common_glued_globalRatio` は
**分母 `s` を 1 つ固定して**「有限個の試験元が同時に大域の比になる」を言う。

★★★★しかし段 E3（チャート写像の全射性）を組むには、
**チャートごとに違う分母 `s_i`** について同時にそれが要る:

* チャート `X_{s_i}` の座標環の生成元 `T_i` を `s_j/s_i` の形にしたい
* ★そのために足す分子の切断は `M^{⊗(n_i+1)}` に住む
* ★★`i` ごとに `n_i` が違うと**族が 1 つの加群に入らない**

★★★★★**本ファイルが単一の `n` を与える**——添字を
「どの分母に属するか」を記録した `d : κ → ι'` つきの `κ` に取り替えるだけでよい。

## ★機構 —— `exists_common_pow` は述語ごとに効く

`exists_common_pow (T : Finset κ) (P : κ → ℕ → Prop) (mono) (ex) : ∃ n, ∀ k ∈ T, P k n`
の `P k n` の中で **`s` を `sfam (d k)` に取り替えるだけ**である。
★単調性（`res_pow_bump`・`hagree_bump`）も存在（`exists_common_exponent`）も
**`k` ごとに局所的**なので、証明は 1 文字も変わらない。

## ★★これで何が繋がるか

★★★分母の族 `{s_i}` と、各 `i` の試験元 `T_i` に対し、
**単一の `n`** と分子 `t_k ∈ Γ(X, sheafify(M^{⊗(n+1)}))` があって

    `g_k = t_k / (s_{d k})^{⊗(n+1)}`（大域の比）

★これで拡大した族 `{(s_i)^{⊗(n+1)}} ∪ {t_k}` が**1 つの加群**に入る。
★★非消失軌跡は `nonVanishing_unit_secPow` で `X_{s_i}` のままなので、
被覆もアフィン性（`§9-915`）も保たれる。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★★★★★★★★★★分母が複数の場合の単一指数 -/

open scoped Classical in
/-- ★★★★★★★★★★**分母が複数でも単一の指数が取れる**。

★`§9-CommonGluedRatio` の `exists_common_exponent_family` を、
「どの分母に属するか」を記録した `d : κ → ι'` つきの添字に一般化したものである。
★★証明は元のものの `s` を `sfam (d k)` に取り替えるだけである
——単調性も存在も `k` ごとに局所的だからである。 -/
theorem exists_common_exponent_family_multi {ι κ ι' : Type} [Fintype ι] (T : Finset κ)
    (M : X.PresheafOfModules) (U : ι → X.Opens)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (sfam : ι' → (M.obj (op ⊤) : Type)) (d : κ → ι')
    (g : ∀ k : κ, (Γ(X, nonVanishing M (sfam (d k))) : Type)) :
    ∃ (n : ℕ) (a : κ → ∀ i, (Γ(X, U i) : Type)), ∀ k ∈ T,
      (∀ i, X.presheaf.map
            (homOfLE (basicOpen_trivValue_le M (U i) (e i) (sfam (d k)))).op (g k)
          * (algebraMap (Γ(X, U i) : Type)
              (Γ(X, X.basicOpen (trivValue M (U i) (e i) (sfam (d k)))) : Type)
              (trivValue M (U i) (e i) (sfam (d k)))) ^ n
        = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) (sfam (d k)))) : Type) (a k i)) ∧
      (∀ i j, X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a k i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a k j)) := by
  have key : ∃ n : ℕ, ∀ k ∈ T, ∃ a : ∀ i, (Γ(X, U i) : Type),
      (∀ i, X.presheaf.map
            (homOfLE (basicOpen_trivValue_le M (U i) (e i) (sfam (d k)))).op (g k)
          * (algebraMap (Γ(X, U i) : Type)
              (Γ(X, X.basicOpen (trivValue M (U i) (e i) (sfam (d k)))) : Type)
              (trivValue M (U i) (e i) (sfam (d k)))) ^ n
        = algebraMap (Γ(X, U i) : Type)
            (Γ(X, X.basicOpen (trivValue M (U i) (e i) (sfam (d k)))) : Type) (a i)) ∧
      (∀ i j, X.presheaf.map (homOfLE (inf_le_left : U i ⊓ U j ≤ U i)).op (a i)
          * (transUnit M (U i ⊓ U j)
              (trivialOfLe M (inf_le_left : U i ⊓ U j ≤ U i) (e i))
              (trivialOfLe M (inf_le_right : U i ⊓ U j ≤ U j) (e j))) ^ n
        = X.presheaf.map (homOfLE (inf_le_right : U i ⊓ U j ≤ U j)).op (a j)) := by
    refine exists_common_pow T _ ?_ ?_
    · rintro k n₁ n₂ hle ⟨a, hden, hag⟩
      obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le hle
      exact ⟨fun i => a i * (trivValue M (U i) (e i) (sfam (d k))) ^ m,
        fun i => res_pow_bump (U i) _ _ (a i) n₁ m (hden i),
        hagree_bump M U e (sfam (d k)) n₁ m a hag⟩
    · intro k _
      obtain ⟨n, a, hden, hag⟩ := exists_common_exponent Finset.univ M U
        (fun i _ => hU i) (fun i _ j _ => hUij i j) e (sfam (d k))
        (fun i => X.presheaf.map
          (homOfLE (basicOpen_trivValue_le M (U i) (e i) (sfam (d k)))).op (g k))
        (fun i _ j _ => by rw [res_trans, res_trans])
      exact ⟨n, a, fun i => hden i (Finset.mem_univ i),
        fun i j => hag i (Finset.mem_univ i) j (Finset.mem_univ j)⟩
  obtain ⟨n, hn⟩ := key
  choose! a ha using hn
  exact ⟨n, a, fun k hk => ha k hk⟩

/-! ## ★★★★★★★★★★★★分母が複数でも同時に大域の比になる -/

open scoped Classical in
/-- ★★★★★★★★★★★★**分母が複数でも同時に大域の比になる** —— 段 E3 の帳簿。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

分母の族 `{s_{i'}}` と、各試験元 `g_k`（分母は `s_{d k}`）に対し、
★**単一の `n`** と分子 `t_k ∈ Γ(X, sheafify(M^{⊗(n+1)}))` があって

    `g_k = t_k / (s_{d k})^{⊗(n+1)}`（大域の比）

★★これで拡大した族が**1 つの加群**に入る。 -/
theorem exists_common_glued_globalRatio_multi {ι κ ι' : Type} [Fintype ι] (T : Finset κ)
    (M : X.PresheafOfModules) (hM : IsLocallyTrivial X M)
    (U : ι → X.Opens) (hcov : (⨆ i, U i) = ⊤)
    (hU : ∀ i, IsAffineOpen (U i)) (hUij : ∀ i j, IsAffineOpen (U i ⊓ U j))
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj M ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (sfam : ι' → (M.obj (op ⊤) : Type)) (d : κ → ι')
    (g : ∀ k : κ, (Γ(X, nonVanishing M (sfam (d k))) : Type)) :
    ∃ (n : ℕ) (t : κ → (((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val.obj
        (op (⊤ : X.Opens)) : Type)), ∀ k ∈ T,
      X.presheaf.map
          (homOfLE (le_of_eq (nonVanishing_unit_secPow M hM (sfam (d k)) n))).op (g k)
        = globalRatio ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val
            (isLocallyTrivial_sheafify X _ (isLocallyTrivial_presheafTensorPow hM (n + 1)))
            (t k) (((sheafifyUnit X (presheafTensorPow M (n + 1))).app (op ⊤)).hom
              (secPow M (sfam (d k)) (n + 1))) := by
  obtain ⟨n, a, hka⟩ := exists_common_exponent_family_multi T M U hU hUij e sfam d g
  have hT : ∀ k ∈ T, ∃ t : (((sheafifyFunctor X).obj
      (presheafTensorPow M (n + 1))).val.obj (op (⊤ : X.Opens)) : Type),
      ∀ i, trivValue ((sheafifyFunctor X).obj (presheafTensorPow M (n + 1))).val (U i)
          (sheafifyTriv (presheafTensorPow M (n + 1)) (tensorPowTriv (e i) (n + 1))) t
        = a k i * (trivValue M (U i) (e i) (sfam (d k))) ^ 1 :=
    fun k hk => exists_glue_bump_trivValue M U hcov e (sfam (d k)) n 1 (a k) (hka k hk).2
  choose! t ht using hT
  exact ⟨n, t, fun k hk =>
    globalRatio_of_den_agree M hM U hcov e (sfam (d k)) n (g k) (a k) (t k)
      (hka k hk).1 (ht k hk)⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_common_exponent_family_multi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分母が複数でも単一の指数が取れる)",
    sectionId := "genell-prop-1-4" }

def exists_common_glued_globalRatio_multi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分母が複数でも同時に大域の比になる)",
    sectionId := "genell-prop-1-4" }

def exists_common_glued_globalRatio_multi.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_common_pow(述語ごとに単一の指数、§9-826)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_common_pow") 1,
    .citation "[ABC3]" "globalRatio_of_den_agree(分子と分母が揃えば大域の比、§9-845)"
      (.inProject "ABC3" "ABC3.Found.GenEll.globalRatio_of_den_agree") 2,
    .implicitStep
      ("★★★★段 E3(チャート写像の全射性)を組むには**チャートごとに違う分母** s_i について" ++
       "同時に指数が要る——i ごとに n_i が違うと族が 1 つの加群に入らないからである。" ++
       "本ファイルが単一の n を与える") 4,
    .implicitStep
      ("★機構: exists_common_pow の述語 P k n の中で s を sfam (d k) に取り替えるだけ。" ++
       "単調性(res_pow_bump・hagree_bump)も存在(exists_common_exponent)も " ++
       "k ごとに局所的なので、証明は 1 文字も変わらない") 2,
    .implicitStep
      ("★★次は拡大した族 {(s_i)^{⊗(n+1)}} ∪ {t_k} を Fin (N+1) に並べ、" ++
       "§9-GlobalChartSurjective の exists_finset_surjective_globalAwayHom に渡す段である。" ++
       "★非消失軌跡は nonVanishing_unit_secPow で X_{s_i} のままなので、" ++
       "被覆もアフィン性(§9-915)も保たれ、埋め込み性は部分族で足りる(§9-916)") 5 ]

end ABC3.Found.GenEll
