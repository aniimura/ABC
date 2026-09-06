import ABC3.Found.PGC.DworkTheta
import Mathlib.RingTheory.PowerSeries.Expand

/-!
# Milne Prop 3.10 の Step 2 —— `σθ ∘ f = g ∘ θ`(`sorry` 無し)

経路 Λ の節点 **Λ6 の最後の段**。`DworkTheta.lean`(Step 1)が与えた
`θ`(`σθ = θ ∘ [u]_f`、`θ(0) = 0`、`θ'(0)` は単元)から出発して、`θ` を取り換え

```
σθ ∘ f = g ∘ θ      (f ∈ F_π, g ∈ F_ϖ, ϖ = uπ)
```

を満たす `θ` を作る(`exists_arithFrobenius_dworkTheta_step2`)。

★**Step 3・Step 4 は Milne THEOREM 3.9 には不要**である。原典自身が
「The proofs of Steps 3 and 4 are straightforward applications of Lemma 2.11,
and so will be left to the reader.」と畳んでおり、THEOREM 3.9 の証明
(「PROOF (THAT K_π · K^ur IS INDEPENDENT OF π)」)が実際に使うのは
`(σθ)(f(T)) = g(θ(T))` だけである。したがって
★**本ファイルが通ったことで Λ6(Dwork による `φ_π` の π 非依存性)の数学は揃う。**

## 原典(Milne, *Class Field Theory*, PROOF (OF STEP 2)、物理 p.51)

```
h = σθ ∘ f ∘ θ⁻¹ = θ ∘ [u]_f ∘ f ∘ θ⁻¹ = θ ∘ f ∘ [u]_f ∘ θ⁻¹.
Then, because f and [u]_f have coefficients in A,
  σh = σθ ∘ f ∘ [u]_f ∘ σ(θ⁻¹) = σθ ∘ f ∘ θ⁻¹ = h.
For the middle equality, we used that [u]_f ∘ σ(θ⁻¹) = θ⁻¹ which follows from
θ ∘ [u]_f ∘ σ(θ⁻¹) = T. Because σh = h, it lies in A[[T]]. Moreover
  h(T) = σε · π · ε⁻¹ T + ⋯ = ϖT + terms of degree ≥ 2,
and h(T) ≡ σθ ∘ (θ⁻¹)^q ≡ (θ(θ⁻¹(T)))^q ≡ T^q mod 𝔪_K. Therefore h ∈ F_ϖ.
Let θ₀ = [1]_{g,h} ∘ θ. Then θ₀ still satisfies (a), and it still satisfies (b)
because [1]_{g,h} ∈ A[[T]]. Moreover σθ₀ ∘ f ∘ θ₀⁻¹ = [1]_{g,h} ∘ h ∘ [1]_{g,h}⁻¹ = g.
```

## 5 つの段と、それぞれの到達点

| 段 | 主張 | 名前 |
|---|---|---|
| 1 | `h := σθ ∘ f ∘ θ'` を作る | (`exists_dworkH_of_residue` の中) |
| 2 | `σh = h` | `map_dworkH_eq_self`(抽象版) |
| 3 | `h ∈ 𝒪_K[[T]]` | `exists_map_eq_of_coeff_mem_range` + 固定環 |
| 4 | `h ∈ F_ϖ`(`h ≡ ϖT (mod deg 2)`・`h ≡ T^q (mod 𝔪)`) | `coeff_one_dworkH` / `map_dworkH_eq_X_pow` |
| 5 | `θ₀ := [1]_{h,g} ∘ θ` で `σθ₀ ∘ f = g ∘ θ₀` | `exists_dworkTheta_step2_of_residue` |

★**段 4 は通った**(落としていない)。★**段 5 も通った**ので、`g` は
「こちらが作った `h`」ではなく**任意に与えられた `F_ϖ` の元**でよい。

## 段 2 で A 有理性が本当に効いたか —— 効いた

`DworkTheta.lean`(Step 1)の逸脱記録 1 が「`f`・`[u]_f` の係数が `𝒪_K` に入る
ことを Step 1 は使っていない」と述べていた。**本ファイルの段 2 はそれを使う。**
実際、抽象版 `map_dworkH_eq_self` の仮定に

* `hFφ : PowerSeries.map φ F = F`
* `hΨφ : PowerSeries.map φ Ψ = Ψ`

の 2 本が入っており、これが「`f` と `[u]_f` の係数が `A` にある」の形式化である。
具体版ではこの 2 本を `map_unramGalCompletionIntHom_map_baseIntHom`
(`σ(ι a) = ι a`)で供給する。**この 2 本を落とすと段 2 は成り立たない**——
`σh = σθ ∘ σf ∘ σ(θ⁻¹)` の中の `σf` を `f` に戻せず、
`[u]_f ∘ σ(θ⁻¹) = θ⁻¹` を使う中央の等式が繋がらない。

## 段 4(`h ≡ T^q mod 𝔪`)の中身

原典の 1 行 `h(T) ≡ σθ ∘ (θ⁻¹)^q ≡ (θ(θ⁻¹(T)))^q ≡ T^q` は、剰余体
`F := 𝓀_{K̂^{ur}}` の上で

```
(P ∘ Z)^q = (map (iterateFrobenius F p f) P) ∘ (Z^q)      (q = p^f)
```

という**標数 `p` の合成版「新入生の夢」**(`pow_subst_eq_subst_map_iterateFrobenius`)
に帰着する。これは mathlib の

* `MvPowerSeries.map_iterateFrobenius_expand`(`map (Frob^n) (expand (p^n) φ) = φ^{p^n}`)
* `PowerSeries.expand_subst`(`(P ∘ Z).expand = P ∘ (Z.expand)`)

の 2 本を継ぐだけで出た。★`FrobeniusExpand.lean` の
`mvPowerSeries_pow_card_eq_expand` は**有限体**専用なのでここでは使えない
(`𝓀_{K̂^{ur}}` は無限)。一般の可換環で成り立つ mathlib 側の補題を直接使う。

段 4 の残りの配管は 3 本。

1. `residue_unramGalCompletionInt`(`DworkAdditive.lean`)——`σ` は
   `𝓀_{K̂^{ur}}` 上で `z ↦ z^q`(Λ5b)。
2. `expChar_residueField_unramifiedCompletionInt`——`𝓀_{K̂^{ur}}` の標数指数が
   `𝓀_K` のそれと同じであること(`IsLocalRing.ResidueField.map (baseIntHom K)` が
   単射なので `charP_of_injective_ringHom` で移す)。
3. `map_residue_map_baseIntHom`——`f ≡ X^q (mod 𝔪_K)` を `𝔪_{K̂^{ur}}` へ移す
   (`isLocalHom_baseIntHom` が既に `DworkFixedRing.lean` にあった)。

## 段 5 の `[1]_{h,g}` の向き

在庫の `LubinTateEndo hq hπmax G hG0 hG1 hG F hF0 hF1 hF a` は関数等式

```
PowerSeries.subst E F = PowerSeries.subst G E      (E := LubinTateEndo …)
```

すなわち `F ∘ E = E ∘ G` を満たす(`LubinTateEndo_functional_equation`)。
段 5 が要るのは `μ ∘ h = g ∘ μ` なので、
★**第 1 引数(`G` の枠)に `h`、第 2 引数(`F` の枠)に `g` を入れる**。
向きを取り違えて `[1]_{g,h}` を作ると `h ∘ μ = μ ∘ g` になり、
`σθ₀ ∘ f = g ∘ θ₀` の計算が閉じない。★Step 1 では `[u]_f` の係数しか
使っていなかったので、**向きの取り違えの初出がここである**。

なお `LubinTateEndo` は 2 本の級数が**同じ素元**を先頭係数に持つことを要求する。
段 5 では `h` も `g` も `F_ϖ`(`ϖ = uπ`)の元なので条件を満たし、
`hϖmax : 𝔪 = (uπ)` は `Ideal.span_singleton_mul_left_unit` で `hπmax` から出る。

## 抽象化の切り分け(`DworkTheta.lean` で効いた手を踏襲した)

段 1–4 の中核は**一般の可換環 `B` と環準同型 `φ : B →+* B`**、および
剰余の側は**一般の可換環 `F` と `r : B →+* F`** の上で書いた。
`PAdicLocalField` の重いインスタンス探索が中核の外に出るので、
抽象版の `lean_check` は **0.05–0.4 秒**で返る。具体版だけが **3–15 秒**である
(`maxHeartbeats 2000000` を付けてある)。★見立てどおり効いた。

## 到達点

| | 主張 |
|---|---|
| 1 | `subst_subst_eq`:合成の結合律(`(P∘a)∘b = P∘(a∘b)`) |
| 2 | `subst_map_substInv_eq`:`Ψ ∘ σθ' = θ'`(原典の「中央の等式」) |
| 3 | `map_dworkH_eq_self`:**`σh = h`**(抽象版、段 2) |
| 4 | `coeff_one_dworkH`:`h'(0) = Ψ'(0)·F'(0)`(段 4 の前半) |
| 5 | `pow_subst_eq_subst_map_iterateFrobenius`:合成版の新入生の夢 |
| 6 | `map_dworkH_eq_X_pow`:**`h ≡ T^q (mod 𝔪)`**(抽象版、段 4 の後半) |
| 7 | `exists_map_eq_of_coeff_mem_range`:係数が像に入る級数は像(段 3) |
| 8 | `map_residue_eq_X_pow_of_map_eq`:剰余体の主張の降下 |
| 9 | `exists_dworkH_of_residue`:**段 1–4 の具体版**(`h ∈ F_ϖ` と `σθ∘f = h∘θ`) |
| 10 | `exists_dworkTheta_step2_of_residue`:**段 5 込みの具体版** |
| 11 | `exists_arithFrobenius_dworkH`:`σ` を `∃` の内側に閉じ込めた段 1–4 |
| 12 | **`exists_arithFrobenius_dworkTheta_step2`**:主定理 |
| 13 | `exists_arithFrobenius_isCoherent_dworkThetaStep2`:Λ5・固定環と同じ 1 つの `σ` |

★**量化の順が主張の中身**である。`∃ σ, ∀ (π,f,u,g), ∃ θ` を主張しており、
`σ` も `θ` も `∃` の内側にある(`∀ (π,f,u,g), ∃ σ` なら `f` ごとに Frobenius を
選び直せるので空虚)。

## 逸脱(記録)

1. **`θ'` を自前で作っている。** 原典は「(a) が `θ` の可逆性を含意する」と書き、
   `DworkTheta.lean` にも `exists_dworkTheta_isIso_of_residue`(合成逆つき)が
   あるが、そちらの結論は `θ'` の**定数項が 0 であること**を出していない。
   段 2 の代入計算には `θ'(0) = 0` が要るので、本ファイルは
   `exists_dworkTheta_of_residue`(合成逆なし)を受けて
   `θ' := θ.substInvOfIsUnit` を自分で置き、mathlib の
   `constantCoeff_substInvOfIsUnit` から `θ'(0) = 0` を取っている。
   ★結論は原典と同じ。段取り係の指示にある「`exists_dworkTheta_isIso_of_residue`
   が入口」から**この 1 点だけずらした**。
2. **`u` の型**(Step 1 と同じ)。原典は `ϖ = uπ` という 2 つの素元の比として
   `u ∈ 𝒪_K^×` を出すが、本ファイルは `u ∈ 𝒪_K` と `IsUnit u` を直接の仮定に
   している。`ϖ` は `u * π` として現れる。★主張は原典より一般。
3. **`g ∈ F_ϖ` を 3 つの等式で書いている**(`coeff 0 g = 0`、
   `coeff 1 g = u * π`、`map residue g = X^q`)。原典の `F_ϖ` の定義そのままで、
   `LubinTateEndo` が要求する形と一致する。
4. **形式群律(条件 (c)・(d))は扱っていない。** 原典 PROPOSITION 3.10 の
   Step 3(`θ(F_f(X,Y)) = F_g(θX,θY)`)と Step 4(`θ ∘ [a]_f = [a]_g ∘ θ`)は
   本ファイルには無い。★冒頭に述べたとおり THEOREM 3.9 には不要である
   (原典自身が「left to the reader」と畳んでいる箇所)。必要になったら
   新しい節点として立てる。
5. **段 2 の抽象版は `Ψ ∘ F = F ∘ Ψ` を仮定の形で受ける**
   (`hcomm`)。具体版では `LubinTateAction_functional_equation` を
   `PowerSeries.map (baseIntHom K)` で送って供給する。

## 退化の自己検査

* **`f`・`[u]_f` の係数が `𝒪_K` に入ることを落とすと段 2 は成り立たない。**
  ★これは Step 1 との違いである(Step 1 は使っていなかった)。上の
  「段 2 で A 有理性が本当に効いたか」を参照。抽象版で言えば `hFφ`・`hΨφ` を
  落とすと `map_dworkH_eq_self` の証明が最初の 2 行で止まる。
* **`σ` を任意の `unramGal K` に替えると偽**。`σ = 1` なら Step 1 の段階で
  `θ = θ ∘ [u]_f` から `u = 1` が出る(`DworkTheta.lean` の自己検査)。
  さらに段 3 も崩れる——`σ = 1` の固定環は `𝒪_{K̂^{ur}}` 全体なので
  「`σh = h` ⇒ `h ∈ 𝒪_K[[T]]`」が言えない。
* **`θ.coeff 1` が単元であることを落とすと `θ'` が作れない**
  (`substInvOfIsUnit` の仮定そのもの)。したがって `h` が定義できない。
* **`u` を非単元まで広げると偽**。`u = 0` なら `[0]_f = 0` で `θ = 0` となり
  (`DworkTheta.lean` の自己検査)、`θ'` が取れない。また `ϖ = 0` は素元でない。
* **`g` の 3 条件のうち `coeff 1 g = u * π` を落とすと段 5 が偽**。
  `LubinTateEndo` は 2 本の級数が同じ素元を先頭係数に持つことを要求する。
  実際 `g'(0) ≠ h'(0)` なら `μ ∘ h = g ∘ μ` を満たす `μ`(`μ'(0) = 1`)は
  次数 1 の比較で `h'(0) = g'(0)` を強いるので存在しない。
* **剰余体の有限性(`Fintype`)は `𝓀_{K̂^{ur}}` には要らない**。要るのは
  `𝓀_K` の側だけで(`Nat.card 𝓀_K = q` と標数指数の移送)、
  `𝓀_{K̂^{ur}}` は無限体のまま扱っている。★`FrobeniusExpand.lean` の
  有限体版を使わなかった理由である。
* `σ` を「位相的生成元」に弱めた場合は**未確認**。この証明は回らない
  (Dwork の補題・固定環がどちらも算術 Frobenius でしか出ていない)が、
  ★**主張自体が偽になるかは確かめていない**。「偽」とは書かない。

## ★見立てとの差(2026-09-06 の実測)

段取り係の見積は **400–800 行**。実測は本ファイル **1005 行**だが、うち
モジュール/節の docstring が **229 行**、各定理の docstring が **171 行**、
空行が 50 行で、**Lean のコード(宣言と証明)は 556 行**である。
主定理の statement 自体が `∃ σ, ∀ (π,f,u,g), ∃ θ θ'` と長い(1 本 30 行)ので、
証明本体はさらに短い。`lean_check` の往復は **32 回**(うち失敗 6 回、
すべて配管:`constantCoeff_map` 不在・`MvPowerSeries` 版と `PowerSeries` 版の
`rw` 不一致 3 回(★`lean-idioms.md` #82 に追加)・`constantCoeff_map_eq_zero`
の一般化漏れ・`set` した冪級数への dot notation(★#83 に追加))。
★**思ったより安い。** 特に

* 段 4(最も込み入ると見られていた)が **mathlib の `expand` 系 2 本**で
  そのまま出た(★`expand` という語を知らないと到達できない命名だった、
  という `FrobeniusExpand.lean` の記録がここで効いた)。
* 段 5 が **在庫の `LubinTateEndo` にそのまま当たった**(向きだけ注意)。

高かったのは具体版の**エラボレーション時間**だけで、
`exists_dworkH_of_residue` 1 本で 12 秒、`exists_dworkTheta_step2_of_residue` で
11 秒かかる(`maxHeartbeats 2000000`)。★配管の値段はここに集中している。

## 設計上の注意(守ったこと)

* **既存の `Found/PGC/*.lean` を書き換えていない**(import のみ)。
* **`Skeleton` / `Interface` を触っていない**。
* **結論に自由なパラメータを出していない**——主定理の型は `K`(と `pp`・`ff`・
  `hq`)にしか依存せず、`σ` も `θ` も `h` も `∃` / 証明の内側にある。
* ★**`frobeniusCompletionInt` は使っていない**(あれは位相的生成元一般)。
  Dwork の補題も固定環も**算術 Frobenius**でしか出ていない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found.GaloisRep
open scoped NNReal Valued

/-! ## 1. 合成の代数 —— 一般の可換環の上で

`PowerSeries.subst a P` は mathlib の約束で `P(a)`、すなわち `P ∘ a` である。
以下ではこの向きを固定して使う。 -/

section Composition

variable {B C : Type*} [CommRing B] [CommRing C]

/-- **合成の結合律**——`(P ∘ a) ∘ b = P ∘ (a ∘ b)`。
mathlib の `substAlgHom_comp_substAlgHom_apply` を `PowerSeries` の言葉に直したもの。 -/
theorem subst_subst_eq {a b f : PowerSeries B} (ha : PowerSeries.HasSubst a)
    (hb : PowerSeries.HasSubst b) :
    PowerSeries.subst b (PowerSeries.subst a f) = PowerSeries.subst (PowerSeries.subst b a) f := by
  have h := PowerSeries.substAlgHom_comp_substAlgHom_apply (R := B) ha hb f
  simpa only [PowerSeries.coe_substAlgHom] using h

/-- 係数写像は定数項 0 を保つ。 -/
theorem constantCoeff_map_eq_zero (φ : B →+* C) {θ : PowerSeries B}
    (hθ0 : PowerSeries.constantCoeff θ = 0) :
    PowerSeries.constantCoeff (PowerSeries.map φ θ) = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, PowerSeries.coeff_map,
    PowerSeries.coeff_zero_eq_constantCoeff_apply, hθ0, map_zero]

/-- `PowerSeries.map_subst` の 1 変数版。★mathlib の版は代入する側を
`MvPowerSeries.map` で書くので、そのままでは `PowerSeries.map` を含む式に
`rw` できない(`Did not find an occurrence of the pattern`)。1 度だけ言い直す。 -/
theorem map_subst_powerSeries (φ : B →+* C) {a f : PowerSeries B}
    (ha : PowerSeries.HasSubst a) :
    PowerSeries.map φ (PowerSeries.subst a f)
      = PowerSeries.subst (PowerSeries.map φ a) (PowerSeries.map φ f) :=
  PowerSeries.map_subst ha f

/-- **合成の 1 次係数**——`(P ∘ a)'(0) = P'(0)·a'(0)`。
`LubinTateUniqueness.coeff_subst_eq_of_order_ge` の `n = 0` の場合。 -/
theorem coeff_one_subst_eq {a f : PowerSeries B} (ha0 : PowerSeries.constantCoeff a = 0)
    (hf0 : PowerSeries.constantCoeff f = 0) :
    PowerSeries.coeff 1 (PowerSeries.subst a f)
      = PowerSeries.coeff 1 f * PowerSeries.coeff 1 a := by
  have hord : ((0 + 1 : ℕ) : ℕ∞) ≤ f.order :=
    PowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hf0
  have h := coeff_subst_eq_of_order_ge (δ := f) (h := a) (π := PowerSeries.coeff 1 a) (n := 0)
    ha0 rfl hord (PowerSeries.HasSubst.of_constantCoeff_zero' ha0)
  simpa using h

/-- `X^k ∘ a = a^k`。 -/
theorem subst_X_pow_eq_pow {a : PowerSeries B} (ha : PowerSeries.HasSubst a) (k : ℕ) :
    PowerSeries.subst a (PowerSeries.X ^ k : PowerSeries B) = a ^ k := by
  have h := map_pow (PowerSeries.substAlgHom (R := B) ha) (PowerSeries.X : PowerSeries B) k
  rw [PowerSeries.coe_substAlgHom] at h
  rw [h, PowerSeries.subst_X ha]

end Composition

/-! ## 2. 段 2 —— `σh = h`(抽象版)

一般の可換環 `B` と環準同型 `φ : B →+* B` の上で、Step 1 の出力
(`θ`・合成逆 `θ'`・`φθ = θ ∘ Ψ`)と `F`・`Ψ` が `φ` で固定されることから
`h := φθ ∘ F ∘ θ'` が `φ` で固定されることを示す。 -/

section Step2Fixed

variable {B : Type*} [CommRing B]

/-- **原典の「中央の等式」**——`Ψ ∘ φθ' = θ'`。

原典 Milne の「For the middle equality, we used that `[u]_f ∘ σ(θ⁻¹) = θ⁻¹`
which follows from `θ ∘ [u]_f ∘ σ(θ⁻¹) = T`」の形式化である。

**証明**。`W := Ψ ∘ φθ'` とおく。`φθ = θ ∘ Ψ` と `φ` が合成と可換なことから
`θ ∘ W = φθ ∘ φθ' = φ(θ ∘ θ') = φ(T) = T`。あとは左から `θ'` を合成して
`W = W ∘ T = W ∘ (θ' ∘ θ) = θ' ∘ (θ ∘ W) = θ' ∘ T = θ'`。 -/
theorem subst_map_substInv_eq (φ : B →+* B) {θ θ' Ψ : PowerSeries B}
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hΨ0 : PowerSeries.constantCoeff Ψ = 0)
    (hinv1 : PowerSeries.subst θ' θ = PowerSeries.X)
    (hinv2 : PowerSeries.subst θ θ' = PowerSeries.X)
    (hθ : PowerSeries.map φ θ = PowerSeries.subst Ψ θ) :
    PowerSeries.subst (PowerSeries.map φ θ') Ψ = θ' := by
  have hθHS : PowerSeries.HasSubst θ := PowerSeries.HasSubst.of_constantCoeff_zero' hθ0
  have hθ'HS : PowerSeries.HasSubst θ' := PowerSeries.HasSubst.of_constantCoeff_zero' hθ'0
  have hΨHS : PowerSeries.HasSubst Ψ := PowerSeries.HasSubst.of_constantCoeff_zero' hΨ0
  have hmθ'0 : PowerSeries.constantCoeff (PowerSeries.map φ θ') = 0 :=
    constantCoeff_map_eq_zero φ hθ'0
  have hmθ'HS : PowerSeries.HasSubst (PowerSeries.map φ θ') :=
    PowerSeries.HasSubst.of_constantCoeff_zero' hmθ'0
  set W := PowerSeries.subst (PowerSeries.map φ θ') Ψ with hWdef
  have hW0 : PowerSeries.constantCoeff W = 0 :=
    PowerSeries.constantCoeff_subst_eq_zero hmθ'0 Ψ hΨ0
  have hWHS : PowerSeries.HasSubst W := PowerSeries.HasSubst.of_constantCoeff_zero' hW0
  have hWθ : PowerSeries.subst W θ = PowerSeries.X := by
    rw [hWdef, ← subst_subst_eq hΨHS hmθ'HS, ← hθ, ← map_subst_powerSeries φ hθ'HS, hinv1,
      PowerSeries.map_X]
  calc W = PowerSeries.subst W PowerSeries.X := (PowerSeries.subst_X hWHS).symm
    _ = PowerSeries.subst W (PowerSeries.subst θ θ') := by rw [hinv2]
    _ = PowerSeries.subst (PowerSeries.subst W θ) θ' := subst_subst_eq hθHS hWHS
    _ = PowerSeries.subst PowerSeries.X θ' := by rw [hWθ]
    _ = θ' := PowerSeries.X_subst θ'

/-- **★★★★★★★★★★★★★★★★(段 2、抽象版)**——`h := φθ ∘ F ∘ θ'` は
`φ` で固定される。

★仮定 `hFφ`・`hΨφ`(`F` と `Ψ` の係数が `φ` の固定部分に入る)が**A 有理性**
そのものであり、**ここで初めて効く**(Step 1 は使っていなかった)。

**証明**。`φ` は係数ごとの作用なので合成と可換する。まず
`φ(φθ) = φ(θ ∘ Ψ) = φθ ∘ φΨ = φθ ∘ Ψ`。したがって

```
φh = φθ' ∘ 側から書くと  φh = (φ(φθ)) ∘ φF ∘ φθ' = (φθ ∘ Ψ) ∘ F ∘ φθ'
   = φθ ∘ (Ψ ∘ F) ∘ φθ' = φθ ∘ (F ∘ Ψ) ∘ φθ'      (hcomm)
   = φθ ∘ F ∘ (Ψ ∘ φθ') = φθ ∘ F ∘ θ' = h.        (subst_map_substInv_eq)
```
-/
theorem map_dworkH_eq_self (φ : B →+* B) {θ θ' F Ψ : PowerSeries B}
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hF0 : PowerSeries.constantCoeff F = 0) (hΨ0 : PowerSeries.constantCoeff Ψ = 0)
    (hinv1 : PowerSeries.subst θ' θ = PowerSeries.X)
    (hinv2 : PowerSeries.subst θ θ' = PowerSeries.X)
    (hFφ : PowerSeries.map φ F = F) (hΨφ : PowerSeries.map φ Ψ = Ψ)
    (hcomm : PowerSeries.subst Ψ F = PowerSeries.subst F Ψ)
    (hθ : PowerSeries.map φ θ = PowerSeries.subst Ψ θ) :
    PowerSeries.map φ (PowerSeries.subst θ' (PowerSeries.subst F (PowerSeries.map φ θ)))
      = PowerSeries.subst θ' (PowerSeries.subst F (PowerSeries.map φ θ)) := by
  have hθ'HS : PowerSeries.HasSubst θ' := PowerSeries.HasSubst.of_constantCoeff_zero' hθ'0
  have hFHS : PowerSeries.HasSubst F := PowerSeries.HasSubst.of_constantCoeff_zero' hF0
  have hΨHS : PowerSeries.HasSubst Ψ := PowerSeries.HasSubst.of_constantCoeff_zero' hΨ0
  have hmθ'0 : PowerSeries.constantCoeff (PowerSeries.map φ θ') = 0 :=
    constantCoeff_map_eq_zero φ hθ'0
  have hmθ'HS : PowerSeries.HasSubst (PowerSeries.map φ θ') :=
    PowerSeries.HasSubst.of_constantCoeff_zero' hmθ'0
  have hkey : PowerSeries.subst (PowerSeries.map φ θ') Ψ = θ' :=
    subst_map_substInv_eq φ hθ0 hθ'0 hΨ0 hinv1 hinv2 hθ
  have hmm : PowerSeries.map φ (PowerSeries.map φ θ)
      = PowerSeries.subst Ψ (PowerSeries.map φ θ) := by
    rw [hθ, map_subst_powerSeries φ hΨHS, hΨφ, ← hθ]
  calc PowerSeries.map φ (PowerSeries.subst θ' (PowerSeries.subst F (PowerSeries.map φ θ)))
      = PowerSeries.subst (PowerSeries.map φ θ')
          (PowerSeries.map φ (PowerSeries.subst F (PowerSeries.map φ θ))) :=
        map_subst_powerSeries φ hθ'HS
    _ = PowerSeries.subst (PowerSeries.map φ θ')
          (PowerSeries.subst F (PowerSeries.map φ (PowerSeries.map φ θ))) := by
        rw [map_subst_powerSeries φ hFHS, hFφ]
    _ = PowerSeries.subst (PowerSeries.map φ θ')
          (PowerSeries.subst F (PowerSeries.subst Ψ (PowerSeries.map φ θ))) := by rw [hmm]
    _ = PowerSeries.subst (PowerSeries.map φ θ')
          (PowerSeries.subst Ψ (PowerSeries.subst F (PowerSeries.map φ θ))) := by
        rw [subst_subst_eq hΨHS hFHS, ← hcomm, ← subst_subst_eq hFHS hΨHS]
    _ = PowerSeries.subst (PowerSeries.subst (PowerSeries.map φ θ') Ψ)
          (PowerSeries.subst F (PowerSeries.map φ θ)) := subst_subst_eq hΨHS hmθ'HS
    _ = PowerSeries.subst θ' (PowerSeries.subst F (PowerSeries.map φ θ)) := by rw [hkey]

/-- `h(0) = 0`。 -/
theorem constantCoeff_dworkH (φ : B →+* B) {θ' F θ : PowerSeries B}
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hF0 : PowerSeries.constantCoeff F = 0) :
    PowerSeries.constantCoeff
        (PowerSeries.subst θ' (PowerSeries.subst F (PowerSeries.map φ θ))) = 0 :=
  PowerSeries.constantCoeff_subst_eq_zero hθ'0 _
    (PowerSeries.constantCoeff_subst_eq_zero hF0 _ (constantCoeff_map_eq_zero φ hθ0))

/-- **段 4 の前半(抽象版)**——`h'(0) = Ψ'(0)·F'(0)`。

原典の `h(T) = σε · π · ε⁻¹ T + ⋯ = ϖT + ⋯` にあたる。実際
`h'(0) = (φθ)'(0)·F'(0)·θ''(0) = (φε)·F'(0)·ε⁻¹` で、`φθ = θ ∘ Ψ` の 1 次係数を
比べると `φε = ε·Ψ'(0)`、さらに `θ ∘ θ' = T` の 1 次係数から `ε·ε⁻¹ = 1`。
★`ε` を名指しせずに `θ.coeff 1` のままで回るので、Step 1 の結論
(`IsUnit (θ.coeff 1)`)だけで足りる。 -/
theorem coeff_one_dworkH (φ : B →+* B) {θ θ' F Ψ : PowerSeries B}
    (hθ0 : PowerSeries.constantCoeff θ = 0) (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hF0 : PowerSeries.constantCoeff F = 0) (hΨ0 : PowerSeries.constantCoeff Ψ = 0)
    (hinv1 : PowerSeries.subst θ' θ = PowerSeries.X)
    (hθ : PowerSeries.map φ θ = PowerSeries.subst Ψ θ) :
    PowerSeries.coeff 1 (PowerSeries.subst θ' (PowerSeries.subst F (PowerSeries.map φ θ)))
      = PowerSeries.coeff 1 Ψ * PowerSeries.coeff 1 F := by
  have hmθ0 : PowerSeries.constantCoeff (PowerSeries.map φ θ) = 0 :=
    constantCoeff_map_eq_zero φ hθ0
  have hinner0 : PowerSeries.constantCoeff (PowerSeries.subst F (PowerSeries.map φ θ)) = 0 :=
    PowerSeries.constantCoeff_subst_eq_zero hF0 _ hmθ0
  have hεinv : PowerSeries.coeff 1 θ * PowerSeries.coeff 1 θ' = 1 := by
    have h := coeff_one_subst_eq hθ'0 hθ0
    rw [hinv1] at h
    simpa using h.symm
  have hε : φ (PowerSeries.coeff 1 θ) = PowerSeries.coeff 1 θ * PowerSeries.coeff 1 Ψ := by
    have h1 := congrArg (PowerSeries.coeff 1) hθ
    rw [PowerSeries.coeff_map, coeff_one_subst_eq hΨ0 hθ0] at h1
    exact h1
  rw [coeff_one_subst_eq hθ'0 hinner0, coeff_one_subst_eq hF0 hmθ0, PowerSeries.coeff_map, hε]
  calc PowerSeries.coeff 1 θ * PowerSeries.coeff 1 Ψ * PowerSeries.coeff 1 F
        * PowerSeries.coeff 1 θ'
      = (PowerSeries.coeff 1 θ * PowerSeries.coeff 1 θ')
        * (PowerSeries.coeff 1 Ψ * PowerSeries.coeff 1 F) := by ring
    _ = PowerSeries.coeff 1 Ψ * PowerSeries.coeff 1 F := by rw [hεinv, one_mul]

end Step2Fixed

/-! ## 3. 段 4 の後半 —— 標数 `p` の合成版「新入生の夢」 -/

section Frobenius

variable {R : Type*} [CommRing R]

/-- **★★★★★★★★合成版の新入生の夢**——標数指数 `p` の可換環の上で

```
(P ∘ Z)^{p^n} = (map (iterateFrobenius R p n) P) ∘ (Z^{p^n}).
```

★これが原典の 1 行 `h(T) ≡ σθ ∘ (θ⁻¹)^q ≡ (θ(θ⁻¹(T)))^q ≡ T^q (mod 𝔪)` の中身
である。証明は mathlib の 2 本を継ぐだけ:

* `PowerSeries.expand_subst`:`(P ∘ Z).expand = P ∘ (Z.expand)`
* `MvPowerSeries.map_iterateFrobenius_expand`:`map (Frob^n) (expand (p^n) W) = W^{p^n}`

★有限体は仮定していない(`FrobeniusExpand.lean` の
`mvPowerSeries_pow_card_eq_expand` は有限体専用なので、`𝓀_{K̂^{ur}}` には
当たらない)。 -/
theorem pow_subst_eq_subst_map_iterateFrobenius {pp : ℕ} (hp : pp ≠ 0) [ExpChar R pp]
    (n : ℕ) {Z P : PowerSeries R} (hZ0 : PowerSeries.constantCoeff Z = 0) :
    (PowerSeries.subst Z P) ^ (pp ^ n)
      = PowerSeries.subst (Z ^ (pp ^ n))
          (PowerSeries.map (iterateFrobenius R pp n) P) := by
  have hZHS : PowerSeries.HasSubst Z := PowerSeries.HasSubst.of_constantCoeff_zero' hZ0
  have hqne : pp ^ n ≠ 0 := pow_ne_zero n hp
  have hexp0 : PowerSeries.constantCoeff (PowerSeries.expand (pp ^ n) hqne Z) = 0 := by
    rw [PowerSeries.constantCoeff_expand, hZ0]
  have hexpHS : PowerSeries.HasSubst (PowerSeries.expand (pp ^ n) hqne Z) :=
    PowerSeries.HasSubst.of_constantCoeff_zero' hexp0
  -- ★`MvPowerSeries.expand` / `MvPowerSeries.map` のままだと `rw` の
  -- パターンが合わないので、`PowerSeries` の言葉で 1 度言い直す(定義は同じ)。
  have h1 : PowerSeries.map (iterateFrobenius R pp n)
      (PowerSeries.expand (pp ^ n) hqne (PowerSeries.subst Z P))
      = (PowerSeries.subst Z P) ^ (pp ^ n) :=
    MvPowerSeries.map_iterateFrobenius_expand (σ := Unit) (R := R) pp hp
      (PowerSeries.subst Z P) n
  have h2 : PowerSeries.map (iterateFrobenius R pp n)
      (PowerSeries.expand (pp ^ n) hqne Z) = Z ^ (pp ^ n) :=
    MvPowerSeries.map_iterateFrobenius_expand (σ := Unit) (R := R) pp hp Z n
  have h3 : PowerSeries.expand (pp ^ n) hqne (PowerSeries.subst Z P)
      = PowerSeries.subst (PowerSeries.expand (pp ^ n) hqne Z) P :=
    PowerSeries.expand_subst (pp ^ n) hqne hZHS P
  rw [← h1, h3, map_subst_powerSeries (iterateFrobenius R pp n) hexpHS, h2]

end Frobenius

section Step2Residue

variable {B F : Type*} [CommRing B] [CommRing F]

/-- `r ∘ φ = (·)^q ∘ r`(Λ5b)を冪級数に持ち上げた形。 -/
theorem map_map_eq_map_iterateFrobenius (r : B →+* F) (φ : B →+* B) {pp ff : ℕ} [ExpChar F pp]
    (hrφ : ∀ b : B, r (φ b) = (r b) ^ (pp ^ ff)) (θ : PowerSeries B) :
    PowerSeries.map r (PowerSeries.map φ θ)
      = PowerSeries.map (iterateFrobenius F pp ff) (PowerSeries.map r θ) := by
  ext n
  rw [PowerSeries.coeff_map, PowerSeries.coeff_map, PowerSeries.coeff_map,
    PowerSeries.coeff_map, hrφ, iterateFrobenius_def]

/-- **★★★★★★★★★★★★(段 4 の後半、抽象版)**——`h ≡ T^q (mod 𝔪)`。

`r : B →+* F` を剰余写像、`q = p^f`、`r(φ b) = (r b)^q`(Λ5b)、
`r F = X^q`(`f ≡ X^q (mod 𝔪_K)`)とすると、`h := φθ ∘ F ∘ θ'` について
`r h = X^q`。

**証明**。`Θ := rθ`、`Y := rθ'` とおくと

```
r h = (r(φθ)) ∘ (rF) ∘ Y = (Frob^f Θ) ∘ X^q ∘ Y = (Frob^f Θ) ∘ Y^q
    = (Θ ∘ Y)^q = (r(θ ∘ θ'))^q = (rT)^q = T^q.
```

3 番目の等号が `pow_subst_eq_subst_map_iterateFrobenius`(合成版の新入生の夢)。
★`θ` の定数項は使っていない(`θ'` の側だけで足りる)。 -/
theorem map_dworkH_eq_X_pow (r : B →+* F) (φ : B →+* B) {pp : ℕ} (hp : pp ≠ 0) [ExpChar F pp]
    {ff : ℕ} (hrφ : ∀ b : B, r (φ b) = (r b) ^ (pp ^ ff))
    {θ θ' Fs : PowerSeries B} (hθ'0 : PowerSeries.constantCoeff θ' = 0)
    (hFs0 : PowerSeries.constantCoeff Fs = 0)
    (hinv1 : PowerSeries.subst θ' θ = PowerSeries.X)
    (hFres : PowerSeries.map r Fs = PowerSeries.X ^ (pp ^ ff)) :
    PowerSeries.map r (PowerSeries.subst θ' (PowerSeries.subst Fs (PowerSeries.map φ θ)))
      = PowerSeries.X ^ (pp ^ ff) := by
  have hθ'HS : PowerSeries.HasSubst θ' := PowerSeries.HasSubst.of_constantCoeff_zero' hθ'0
  have hFsHS : PowerSeries.HasSubst Fs := PowerSeries.HasSubst.of_constantCoeff_zero' hFs0
  set Y := PowerSeries.map r θ' with hYdef
  set Θ := PowerSeries.map r θ with hΘdef
  have hY0 : PowerSeries.constantCoeff Y = 0 := constantCoeff_map_eq_zero r hθ'0
  have hYHS : PowerSeries.HasSubst Y := PowerSeries.HasSubst.of_constantCoeff_zero' hY0
  have hXqHS : PowerSeries.HasSubst ((PowerSeries.X : PowerSeries F) ^ (pp ^ ff)) := by
    apply PowerSeries.HasSubst.of_constantCoeff_zero'
    rw [map_pow, PowerSeries.constantCoeff_X, zero_pow (pow_ne_zero ff hp)]
  have hYΘ : PowerSeries.subst Y Θ = PowerSeries.X := by
    rw [hYdef, hΘdef, ← map_subst_powerSeries r hθ'HS, hinv1, PowerSeries.map_X]
  calc PowerSeries.map r (PowerSeries.subst θ' (PowerSeries.subst Fs (PowerSeries.map φ θ)))
      = PowerSeries.subst Y
          (PowerSeries.map r (PowerSeries.subst Fs (PowerSeries.map φ θ))) :=
        map_subst_powerSeries r hθ'HS
    _ = PowerSeries.subst Y (PowerSeries.subst ((PowerSeries.X : PowerSeries F) ^ (pp ^ ff))
          (PowerSeries.map (iterateFrobenius F pp ff) Θ)) := by
        rw [map_subst_powerSeries r hFsHS, hFres, map_map_eq_map_iterateFrobenius r φ hrφ θ]
    _ = PowerSeries.subst
          (PowerSeries.subst Y ((PowerSeries.X : PowerSeries F) ^ (pp ^ ff)))
          (PowerSeries.map (iterateFrobenius F pp ff) Θ) := subst_subst_eq hXqHS hYHS
    _ = PowerSeries.subst (Y ^ (pp ^ ff)) (PowerSeries.map (iterateFrobenius F pp ff) Θ) := by
        rw [subst_X_pow_eq_pow hYHS]
    _ = (PowerSeries.subst Y Θ) ^ (pp ^ ff) :=
        (pow_subst_eq_subst_map_iterateFrobenius hp ff hY0).symm
    _ = PowerSeries.X ^ (pp ^ ff) := by rw [hYΘ]

end Step2Residue

/-! ## 4. 段 3 —— `h ∈ ι(A[[T]])` への降下 -/

section Descent

variable {A B : Type*} [CommRing A] [CommRing B]

/-- 係数がすべて `ι` の像に入る冪級数は `ι` の像である。★段 3 の中身は
これだけで、`𝔪` 進完備性も収束も使わない(係数ごとに選ぶ)。 -/
theorem exists_map_eq_of_coeff_mem_range (ι : A →+* B) {h : PowerSeries B}
    (hc : ∀ n : ℕ, ∃ a : A, ι a = PowerSeries.coeff n h) :
    ∃ H : PowerSeries A, PowerSeries.map ι H = h := by
  classical
  choose g hg using hc
  refine ⟨PowerSeries.mk g, ?_⟩
  ext n
  rw [PowerSeries.coeff_map, PowerSeries.coeff_mk, hg]

/-- 剰余体での主張の降下——`ι` が局所準同型なら、剰余体の写像
`IsLocalRing.ResidueField.map ι` は体からの環準同型なので単射であり、
`h = ι(H)` の剰余の情報は `H` に降りる。 -/
theorem map_residue_eq_X_pow_of_map_eq [IsLocalRing A] [IsLocalRing B]
    (ι : A →+* B) [IsLocalHom ι] {H : PowerSeries A} {h : PowerSeries B}
    (hH : PowerSeries.map ι H = h) {q : ℕ}
    (hres : PowerSeries.map (IsLocalRing.residue B) h = PowerSeries.X ^ q) :
    PowerSeries.map (IsLocalRing.residue A) H = PowerSeries.X ^ q := by
  set j := IsLocalRing.ResidueField.map ι with hjdef
  have hjinj : Function.Injective j := j.injective
  apply PowerSeries.map_injective j hjinj
  have hcomm : PowerSeries.map j (PowerSeries.map (IsLocalRing.residue A) H)
      = PowerSeries.map (IsLocalRing.residue B) (PowerSeries.map ι H) := by
    ext n
    rw [PowerSeries.coeff_map, PowerSeries.coeff_map, PowerSeries.coeff_map,
      PowerSeries.coeff_map, hjdef, IsLocalRing.ResidueField.map_residue]
  rw [hcomm, hH, hres, map_pow, PowerSeries.map_X]

end Descent

/-! ## 5. 具体版 —— `B = 𝒪_{K̂^{ur}}`、`φ = σ`(算術 Frobenius) -/

variable {p : ℕ} [Fact p.Prime]

/-- 有限かつ非自明な環では、標数指数は素数であり、実際に標数である。
(`ExpChar` のもう一方の枝は `CharZero` を含み、有限性と両立しない。) -/
theorem prime_of_expChar_of_fintype (R : Type*) [CommRing R] [Fintype R] [Nontrivial R]
    (pp : ℕ) [ExpChar R pp] : pp.Prime ∧ CharP R pp := by
  rcases ‹ExpChar R pp› with _ | ⟨hprime⟩
  · have : Infinite R := inferInstance
    exact absurd (not_finite R) (by simp)
  · exact ⟨hprime, ‹_›⟩

/-- `𝓀_{K̂^{ur}}` の標数指数は `𝓀_K` のそれと同じ——`𝓀_K → 𝓀_{K̂^{ur}}` が
単射だから(`charP_of_injective_ringHom`)。★`𝓀_{K̂^{ur}}` は無限体でよい。 -/
theorem expChar_residueField_unramifiedCompletionInt (K : PAdicLocalField p) {pp : ℕ}
    [ExpChar 𝓀[K.carrier] pp] [Fintype 𝓀[K.carrier]] :
    ExpChar (IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)) pp := by
  obtain ⟨hprime, hchar⟩ := prime_of_expChar_of_fintype 𝓀[K.carrier] pp
  haveI := hchar
  haveI : CharP (IsLocalRing.ResidueField ↥(unramifiedCompletionInt K)) pp :=
    charP_of_injective_ringHom
      (IsLocalRing.ResidueField.map (baseIntHom K)).injective pp
  exact ExpChar.prime hprime

/-- **★A 有理性**——`𝒪_K` 係数の冪級数は `σ` で動かない。
段 2 の `hFφ`・`hΨφ` を供給するのがこの補題である。 -/
theorem map_unramGalCompletionIntHom_map_baseIntHom (K : PAdicLocalField p) (σ : unramGal K)
    (a : PowerSeries 𝒪[K.carrier]) :
    PowerSeries.map (unramGalCompletionIntHom K σ) (PowerSeries.map (baseIntHom K) a)
      = PowerSeries.map (baseIntHom K) a := by
  ext n
  rw [PowerSeries.coeff_map, PowerSeries.coeff_map, unramGalCompletionIntHom_apply,
    unramGalCompletionInt_baseIntHom]

/-- `f ≡ X^q (mod 𝔪_K)` を `𝔪_{K̂^{ur}}` の側へ移す。 -/
theorem map_residue_map_baseIntHom (K : PAdicLocalField p) {q : ℕ}
    {f : PowerSeries 𝒪[K.carrier]}
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ q) :
    PowerSeries.map (IsLocalRing.residue ↥(unramifiedCompletionInt K))
        (PowerSeries.map (baseIntHom K) f) = PowerSeries.X ^ q := by
  have hcomm : PowerSeries.map (IsLocalRing.residue ↥(unramifiedCompletionInt K))
      (PowerSeries.map (baseIntHom K) f)
      = PowerSeries.map (IsLocalRing.ResidueField.map (baseIntHom K))
          (PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f) := by
    ext n
    rw [PowerSeries.coeff_map, PowerSeries.coeff_map, PowerSeries.coeff_map,
      PowerSeries.coeff_map, IsLocalRing.ResidueField.map_residue]
  rw [hcomm, hf, map_pow, PowerSeries.map_X]

set_option maxHeartbeats 2000000 in
/-- **★★★★★★★★★★★★★★★★★★★★(Λ6、Milne Prop 3.10 Step 2 の段 1–4)**——
剰余体上で `q` 乗として作用する `σ ∈ Gal(K^{ur}/K)` と、素元 `π`・
Lubin-Tate 級数 `f ∈ F_π`・単元 `u ∈ 𝒪_K^×` に対し、

```
∃ θ θ' ∈ 𝒪_{K̂^{ur}}[[T]], ∃ h ∈ 𝒪_K[[T]],
  θ(0) = 0 ∧ θ'(0) は単元 ∧ θ ∘ θ' = θ' ∘ θ = T ∧ σθ = θ ∘ [u]_f ∧
  h ∈ F_{uπ} ∧ σθ ∘ f = h ∘ θ.
```

★`h ∈ F_{uπ}` は 3 本の等式(`h(0) = 0`、`h'(0) = uπ`、`h ≡ T^q (mod 𝔪)`)で
書いてある。★**段 4 まで通っている**(落としていない)。

**証明の骨格**(原典 Milne の Step 2 そのまま)。

1. Step 1(`exists_dworkTheta_of_residue`)で `θ` を取り、
   `θ' := θ.substInvOfIsUnit` を置く(逸脱 1)。
2. `h := σθ ∘ f ∘ θ'` に `map_dworkH_eq_self` を当てて `σh = h`。
   ★ここで `f`・`[u]_f` の A 有理性(`map_unramGalCompletionIntHom_map_baseIntHom`)と
   関数等式 `f ∘ [u]_f = [u]_f ∘ f`(`LubinTateAction_functional_equation`)を使う。
3. 固定環(`unramGalCompletionInt_eq_self_iff`)を係数ごとに当て、
   `exists_map_eq_of_coeff_mem_range` で `h = ι(H)` と書く。
4. `coeff_one_dworkH` で `H'(0) = uπ`、`map_dworkH_eq_X_pow` +
   `map_residue_eq_X_pow_of_map_eq` で `H ≡ T^q (mod 𝔪_K)`。
5. 最後の等式は `h ∘ θ = (σθ ∘ f ∘ θ') ∘ θ = σθ ∘ f`(結合律と `θ' ∘ θ = T`)。 -/
theorem exists_dworkH_of_residue (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u) :
    ∃ (θ θ' : PowerSeries ↥(unramifiedCompletionInt K)) (H : PowerSeries 𝒪[K.carrier]),
      PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
      PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X ∧
      PowerSeries.map (unramGalCompletionIntHom K σ) θ
        = PowerSeries.subst (PowerSeries.map (baseIntHom K)
            (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ ∧
      PowerSeries.coeff 0 H = 0 ∧ PowerSeries.coeff 1 H = u * π ∧
      PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) H = PowerSeries.X ^ (pp ^ ff) ∧
      PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
          (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
        = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) H) := by
  classical
  obtain ⟨hprime, -⟩ := prime_of_expChar_of_fintype (IsLocalRing.ResidueField (𝒪[K.carrier])) pp
  haveI := expChar_residueField_unramifiedCompletionInt (K := K) (pp := pp)
  obtain ⟨θ, h0, h1, heq⟩ := exists_dworkTheta_of_residue K hσ hq hπmax f hf0 hf1 hf hu
  have hθ0 : PowerSeries.constantCoeff θ = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact h0
  refine ⟨θ, θ.substInvOfIsUnit h1, ?_⟩
  set θ' := θ.substInvOfIsUnit h1 with hθ'def
  have hθ'0 : PowerSeries.constantCoeff θ' = 0 :=
    PowerSeries.constantCoeff_substInvOfIsUnit θ h1
  have hinv1 : PowerSeries.subst θ' θ = PowerSeries.X :=
    PowerSeries.subst_substInvOfIsUnit_right θ hθ0 h1
  have hinv2 : PowerSeries.subst θ θ' = PowerSeries.X :=
    PowerSeries.subst_substInvOfIsUnit_left θ hθ0 h1
  have hθHS : PowerSeries.HasSubst θ := PowerSeries.HasSubst.of_constantCoeff_zero' hθ0
  have hθ'HS : PowerSeries.HasSubst θ' := PowerSeries.HasSubst.of_constantCoeff_zero' hθ'0
  set Λ := LubinTateAction hq hπmax f hf0 hf1 hf u with hΛdef
  set Fs := PowerSeries.map (baseIntHom K) f with hFsdef
  set Ψ := PowerSeries.map (baseIntHom K) Λ with hΨdef
  have hf0' : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hΛ0 : PowerSeries.constantCoeff Λ = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf u
  have hFs0 : PowerSeries.constantCoeff Fs = 0 := constantCoeff_map_eq_zero _ hf0'
  have hΨ0 : PowerSeries.constantCoeff Ψ = 0 := constantCoeff_map_eq_zero _ hΛ0
  have hFsφ : PowerSeries.map (unramGalCompletionIntHom K σ) Fs = Fs :=
    map_unramGalCompletionIntHom_map_baseIntHom K σ f
  have hΨφ : PowerSeries.map (unramGalCompletionIntHom K σ) Ψ = Ψ :=
    map_unramGalCompletionIntHom_map_baseIntHom K σ Λ
  have hcomm : PowerSeries.subst Ψ Fs = PowerSeries.subst Fs Ψ := by
    rw [hΨdef, hFsdef,
      ← map_subst_powerSeries (baseIntHom K) (PowerSeries.HasSubst.of_constantCoeff_zero' hΛ0),
      LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf u,
      map_subst_powerSeries (baseIntHom K) (PowerSeries.HasSubst.of_constantCoeff_zero' hf0')]
  have hfix := map_dworkH_eq_self (unramGalCompletionIntHom K σ) hθ0 hθ'0 hFs0 hΨ0 hinv1 hinv2
    hFsφ hΨφ hcomm heq
  set hser := PowerSeries.subst θ'
    (PowerSeries.subst Fs (PowerSeries.map (unramGalCompletionIntHom K σ) θ)) with hserdef
  have hcoeffs : ∀ n : ℕ, ∃ a : 𝒪[K.carrier], baseIntHom K a = PowerSeries.coeff n hser := by
    intro n
    refine (unramGalCompletionInt_eq_self_iff K hσ _).mp ?_
    have hc := congrArg (PowerSeries.coeff n) hfix
    rwa [PowerSeries.coeff_map, unramGalCompletionIntHom_apply] at hc
  obtain ⟨H, hH⟩ := exists_map_eq_of_coeff_mem_range (baseIntHom K) hcoeffs
  have hH0 : PowerSeries.coeff 0 H = 0 := by
    apply baseIntHom_injective K
    rw [← PowerSeries.coeff_map, hH, map_zero, hserdef,
      PowerSeries.coeff_zero_eq_constantCoeff_apply]
    exact constantCoeff_dworkH (unramGalCompletionIntHom K σ) hθ0 hθ'0 hFs0
  have hH1 : PowerSeries.coeff 1 H = u * π := by
    apply baseIntHom_injective K
    rw [← PowerSeries.coeff_map, hH, map_mul, hserdef,
      coeff_one_dworkH (unramGalCompletionIntHom K σ) hθ0 hθ'0 hFs0 hΨ0 hinv1 heq,
      hΨdef, hFsdef, PowerSeries.coeff_map, PowerSeries.coeff_map, hΛdef,
      coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf u, hf1]
  have hHres : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) H
      = PowerSeries.X ^ (pp ^ ff) := by
    refine map_residue_eq_X_pow_of_map_eq (baseIntHom K) hH ?_
    rw [hserdef]
    refine map_dworkH_eq_X_pow (IsLocalRing.residue _) (unramGalCompletionIntHom K σ)
      hprime.ne_zero ?_ hθ'0 hFs0 hinv1 ?_
    · intro b
      rw [unramGalCompletionIntHom_apply, residue_unramGalCompletionInt K hσ b,
        Nat.card_eq_fintype_card, hq]
    · rw [hFsdef]; exact map_residue_map_baseIntHom K hf
  have hfinal : PowerSeries.subst Fs (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
      = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) H) := by
    rw [hH, hserdef, subst_subst_eq hθ'HS hθHS, hinv2, PowerSeries.X_subst]
  exact ⟨H, h0, h1, hinv1, hinv2, heq, hH0, hH1, hHres, hfinal⟩

set_option maxHeartbeats 2000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★(Λ6、Milne Prop 3.10 Step 2 の全体)**——
`g` を**任意に与えられた** `F_ϖ`(`ϖ = uπ`)の元とするとき、Step 1 の `θ` を
取り換えて `σθ ∘ f = g ∘ θ` にできる。

**証明**(原典の段 5)。`exists_dworkH_of_residue` の `h` と与えられた `g` は
どちらも `F_ϖ` の元なので、在庫の `LubinTateEndo` が `μ := [1]_{h,g}`
(`μ(0) = 0`、`μ'(0) = 1`、`μ ∘ h = g ∘ μ`)を与える。
★**引数の向き**:`LubinTateEndo hq hϖmax h … g … 1` の第 1 引数(`G` の枠)が
`h`、第 2 引数(`F` の枠)が `g` である(モジュール docstring 参照)。
`θ₀ := μ ∘ θ` とおくと

* `θ₀(0) = 0`、`θ₀'(0) = μ'(0)·θ'(0) = θ'(0)` は単元;
* `σθ₀ = σμ ∘ σθ = μ ∘ θ ∘ [u]_f = θ₀ ∘ [u]_f`(`μ` の係数は `𝒪_K`);
* `σθ₀ ∘ f = μ ∘ (σθ ∘ f) = μ ∘ h ∘ θ = g ∘ μ ∘ θ = g ∘ θ₀`。

★`hg1 : g'(0) = u * π` を落とすと `LubinTateEndo` が使えない(2 本の級数が
同じ素元を先頭係数に持つことを要求する)。 -/
theorem exists_dworkTheta_step2_of_residue (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    {u : 𝒪[K.carrier]} (hu : IsUnit u)
    (g : PowerSeries (𝒪[K.carrier])) (hg0 : PowerSeries.coeff 0 g = 0)
    (hg1 : PowerSeries.coeff 1 g = u * π)
    (hg : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff)) :
    ∃ θ θ' : PowerSeries ↥(unramifiedCompletionInt K),
      PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
      PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X ∧
      PowerSeries.map (unramGalCompletionIntHom K σ) θ
        = PowerSeries.subst (PowerSeries.map (baseIntHom K)
            (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ ∧
      PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
          (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
        = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g) := by
  classical
  obtain ⟨θ, θ', H, h0, h1, hinv1, hinv2, heq, hH0, hH1, hHres, hfinal⟩ :=
    exists_dworkH_of_residue K hσ hq hπmax f hf0 hf1 hf hu
  have hϖmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {u * π} := by
    rw [hπmax, Ideal.span_singleton_mul_left_unit hu π]
  have hH0c : PowerSeries.constantCoeff H = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hH0
  set μ := LubinTateEndo hq hϖmax H hH0c hH1 hHres g hg0 hg1 hg 1 with hμdef
  have hμ0 : PowerSeries.constantCoeff μ = 0 :=
    constantCoeff_LubinTateEndo hq hϖmax H hH0c hH1 hHres g hg0 hg1 hg 1
  have hμ1 : PowerSeries.coeff 1 μ = 1 :=
    coeff_one_LubinTateEndo hq hϖmax H hH0c hH1 hHres g hg0 hg1 hg 1
  have hμfe : PowerSeries.subst μ g = PowerSeries.subst H μ :=
    LubinTateEndo_functional_equation hq hϖmax H hH0c hH1 hHres g hg0 hg1 hg 1
  set Fs := PowerSeries.map (baseIntHom K) f with hFsdef
  set Ψ := PowerSeries.map (baseIntHom K) (LubinTateAction hq hπmax f hf0 hf1 hf u) with hΨdef
  set M := PowerSeries.map (baseIntHom K) μ with hMdef
  have hθ0 : PowerSeries.constantCoeff θ = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact h0
  have hθHS : PowerSeries.HasSubst θ := PowerSeries.HasSubst.of_constantCoeff_zero' hθ0
  have hf0' : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hFs0 : PowerSeries.constantCoeff Fs = 0 := constantCoeff_map_eq_zero _ hf0'
  have hFsHS : PowerSeries.HasSubst Fs := PowerSeries.HasSubst.of_constantCoeff_zero' hFs0
  have hΨ0 : PowerSeries.constantCoeff Ψ = 0 :=
    constantCoeff_map_eq_zero _ (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf u)
  have hΨHS : PowerSeries.HasSubst Ψ := PowerSeries.HasSubst.of_constantCoeff_zero' hΨ0
  have hM0 : PowerSeries.constantCoeff M = 0 := constantCoeff_map_eq_zero _ hμ0
  have hMHS : PowerSeries.HasSubst M := PowerSeries.HasSubst.of_constantCoeff_zero' hM0
  have hM1 : PowerSeries.coeff 1 M = 1 := by rw [hMdef, PowerSeries.coeff_map, hμ1, map_one]
  have hMφ : PowerSeries.map (unramGalCompletionIntHom K σ) M = M :=
    map_unramGalCompletionIntHom_map_baseIntHom K σ μ
  have hHs0 : PowerSeries.constantCoeff (PowerSeries.map (baseIntHom K) H) = 0 :=
    constantCoeff_map_eq_zero _ hH0c
  have hHsHS : PowerSeries.HasSubst (PowerSeries.map (baseIntHom K) H) :=
    PowerSeries.HasSubst.of_constantCoeff_zero' hHs0
  have hμHS : PowerSeries.HasSubst μ := PowerSeries.HasSubst.of_constantCoeff_zero' hμ0
  have hHHS : PowerSeries.HasSubst H := PowerSeries.HasSubst.of_constantCoeff_zero' hH0c
  set θ0 : PowerSeries ↥(unramifiedCompletionInt K) := PowerSeries.subst θ M with hθ0def
  have hθ00 : PowerSeries.constantCoeff θ0 = 0 :=
    PowerSeries.constantCoeff_subst_eq_zero hθ0 M hM0
  have hθ01 : PowerSeries.coeff 1 θ0 = PowerSeries.coeff 1 θ := by
    rw [hθ0def, coeff_one_subst_eq hθ0 hM0, hM1, one_mul]
  have hθ01u : IsUnit (PowerSeries.coeff 1 θ0) := by rw [hθ01]; exact h1
  have hmapθ0 : PowerSeries.map (unramGalCompletionIntHom K σ) θ0
      = PowerSeries.subst (PowerSeries.map (unramGalCompletionIntHom K σ) θ) M := by
    rw [hθ0def, map_subst_powerSeries (unramGalCompletionIntHom K σ) hθHS, hMφ]
  refine ⟨θ0, PowerSeries.substInvOfIsUnit θ0 hθ01u, ?_, hθ01u,
    PowerSeries.subst_substInvOfIsUnit_right θ0 hθ00 hθ01u,
    PowerSeries.subst_substInvOfIsUnit_left θ0 hθ00 hθ01u, ?_, ?_⟩
  · rw [PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hθ00
  · calc PowerSeries.map (unramGalCompletionIntHom K σ) θ0
        = PowerSeries.subst (PowerSeries.subst Ψ θ) M := by rw [hmapθ0, heq]
      _ = PowerSeries.subst Ψ θ0 := by rw [hθ0def, subst_subst_eq hθHS hΨHS]
  · calc PowerSeries.subst Fs (PowerSeries.map (unramGalCompletionIntHom K σ) θ0)
        = PowerSeries.subst Fs
            (PowerSeries.subst (PowerSeries.map (unramGalCompletionIntHom K σ) θ) M) := by
          rw [hmapθ0]
      _ = PowerSeries.subst
            (PowerSeries.subst Fs (PowerSeries.map (unramGalCompletionIntHom K σ) θ)) M :=
          subst_subst_eq (PowerSeries.HasSubst.of_constantCoeff_zero'
            (constantCoeff_map_eq_zero _ hθ0)) hFsHS
      _ = PowerSeries.subst (PowerSeries.subst θ (PowerSeries.map (baseIntHom K) H)) M := by
          rw [hfinal]
      _ = PowerSeries.subst θ (PowerSeries.subst (PowerSeries.map (baseIntHom K) H) M) :=
          (subst_subst_eq hHsHS hθHS).symm
      _ = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) (PowerSeries.subst H μ)) := by
          rw [map_subst_powerSeries (baseIntHom K) hHHS]
      _ = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) (PowerSeries.subst μ g)) := by
          rw [hμfe]
      _ = PowerSeries.subst θ (PowerSeries.subst M (PowerSeries.map (baseIntHom K) g)) := by
          rw [map_subst_powerSeries (baseIntHom K) hμHS]
      _ = PowerSeries.subst θ0 (PowerSeries.map (baseIntHom K) g) := by
          rw [hθ0def, subst_subst_eq hMHS hθHS]

/-! ## 6. 主定理 —— `σ` を `∃` の内側に閉じ込めた形 -/

set_option maxHeartbeats 2000000 in
/-- **★★★★★★★★★★★★★★★★★★★★(Λ6、段 1–4 の主定理)**——
`Gal(K^{ur}/K)` には、剰余体上で `z ↦ z^q` として作用し(= 算術 Frobenius)、
**その 1 つの `σ` について**、任意の `(π, f, u)` に対し
`h ∈ F_{uπ} ⊆ 𝒪_K[[T]]` と `θ` が取れて `σθ ∘ f = h ∘ θ` が成り立つ元がある。 -/
theorem exists_arithFrobenius_dworkH (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      ∀ (π : 𝒪[K.carrier])
        (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
        (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
        (hf1 : PowerSeries.coeff 1 f = π)
        (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
        (u : 𝒪[K.carrier]), IsUnit u →
        ∃ (θ θ' : PowerSeries ↥(unramifiedCompletionInt K)) (H : PowerSeries 𝒪[K.carrier]),
          PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
          PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X ∧
          PowerSeries.map (unramGalCompletionIntHom K σ) θ
            = PowerSeries.subst (PowerSeries.map (baseIntHom K)
                (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ ∧
          PowerSeries.coeff 0 H = 0 ∧ PowerSeries.coeff 1 H = u * π ∧
          PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) H = PowerSeries.X ^ (pp ^ ff) ∧
          PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
              (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
            = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) H) :=
  ⟨arithFrobenius K, arithFrobenius_residue K,
    fun _π hπmax f hf0 hf1 hf _u hu =>
      exists_dworkH_of_residue K (arithFrobenius_residue K) hq hπmax f hf0 hf1 hf hu⟩

set_option maxHeartbeats 2000000 in
/-- **★★★★★★★★★★★★★★★★★★★★★★★★★★(Λ6、Milne Prop 3.10 Step 2 の主定理)**——
`Gal(K^{ur}/K)` には、

* `𝓀_{K^{ur}}` 上で `z ↦ z^q` として作用し(= 算術 Frobenius)、
* **その 1 つの `σ` について**、任意の素元 `π`・任意の `f ∈ F_π`・
  任意の単元 `u ∈ 𝒪_K^×`・任意の `g ∈ F_{uπ}` に対し、
  `θ(0) = 0`、`θ'(0)` が単元、`θ` は合成に関して可逆、`σθ = θ ∘ [u]_f`、
  そして
  ```
  σθ ∘ f = g ∘ θ
  ```
  なる `θ ∈ 𝒪_{K̂^{ur}}[[T]]` が存在する

元 `σ` がある。

★これが Milne THEOREM 3.9(`K_π · K^{ur}` と `φ_π` の `π` 非依存性)が
消費する形である。THEOREM 3.9 の証明は本定理から

```
(σθ)([π]_f(T)) = (σθ)(f(T))     (f = [π]_f)
               = g(θ(T)) = [ϖ]_g(θ(T))
```

を得て、`f` の根と `g` の根を対応させる。★Step 3・Step 4(形式群律との
両立)は THEOREM 3.9 には要らない。

**証明**:`arithFrobenius K` に `exists_dworkTheta_step2_of_residue` を当てる。

退化の自己検査(詳細はモジュール docstring)。

* **`f`・`[u]_f` の係数が `𝒪_K` に入ることを落とすと段 2 が偽**
  (★Step 1 との違い)。
* **`σ` を任意の `unramGal K` に替えると偽**(`σ = 1` なら Step 1 で `u = 1`、
  段 3 でも固定環が全体になって降下できない)。
* **`θ.coeff 1` が単元であることを落とすと `θ'` が作れない**。
* **`hg1`(`g'(0) = uπ`)を落とすと段 5 が偽**——`[1]_{h,g}` が存在しない。 -/
theorem exists_arithFrobenius_dworkTheta_step2 (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      ∀ (π : 𝒪[K.carrier])
        (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
        (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
        (hf1 : PowerSeries.coeff 1 f = π)
        (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
        (u : 𝒪[K.carrier]), IsUnit u →
        ∀ (g : PowerSeries (𝒪[K.carrier])), PowerSeries.coeff 0 g = 0 →
        PowerSeries.coeff 1 g = u * π →
        PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff) →
        ∃ θ θ' : PowerSeries ↥(unramifiedCompletionInt K),
          PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
          PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X ∧
          PowerSeries.map (unramGalCompletionIntHom K σ) θ
            = PowerSeries.subst (PowerSeries.map (baseIntHom K)
                (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ ∧
          PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
              (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
            = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g) :=
  ⟨arithFrobenius K, arithFrobenius_residue K,
    fun _π hπmax f hf0 hf1 hf _u hu g hg0 hg1 hg =>
      exists_dworkTheta_step2_of_residue K (arithFrobenius_residue K) hq hπmax f hf0 hf1 hf hu
        g hg0 hg1 hg⟩

set_option maxHeartbeats 2000000 in
/-- **★★★★★★★★★★★★★★★★★★★★同じ `σ` が Λ5 の位相的生成元でも
固定環の Frobenius でも Step 2 の解でもある形**。

Λ5(`Ẑ` との同定)・Λ6 の固定環・そして本ファイルの Step 2 が**同じ 1 つの元**を
指せるように、4 つの性質を 1 つの `∃` に並べた形。★4 つを別々の `∃` に
分けたら無意味になる。 -/
theorem exists_arithFrobenius_isCoherent_dworkThetaStep2 (K : PAdicLocalField p)
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) ∧
      (∀ b : ↥(unramifiedCompletionInt K),
        unramGalCompletionInt K σ b = b ↔ ∃ a : 𝒪[K.carrier], baseIntHom K a = b) ∧
      ∀ (π : 𝒪[K.carrier])
        (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
        (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
        (hf1 : PowerSeries.coeff 1 f = π)
        (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
        (u : 𝒪[K.carrier]), IsUnit u →
        ∀ (g : PowerSeries (𝒪[K.carrier])), PowerSeries.coeff 0 g = 0 →
        PowerSeries.coeff 1 g = u * π →
        PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) g = PowerSeries.X ^ (pp ^ ff) →
        ∃ θ θ' : PowerSeries ↥(unramifiedCompletionInt K),
          PowerSeries.coeff 0 θ = 0 ∧ IsUnit (PowerSeries.coeff 1 θ) ∧
          PowerSeries.subst θ' θ = PowerSeries.X ∧ PowerSeries.subst θ θ' = PowerSeries.X ∧
          PowerSeries.map (unramGalCompletionIntHom K σ) θ
            = PowerSeries.subst (PowerSeries.map (baseIntHom K)
                (LubinTateAction hq hπmax f hf0 hf1 hf u)) θ ∧
          PowerSeries.subst (PowerSeries.map (baseIntHom K) f)
              (PowerSeries.map (unramGalCompletionIntHom K σ) θ)
            = PowerSeries.subst θ (PowerSeries.map (baseIntHom K) g) :=
  ⟨arithFrobenius K, arithFrobenius_residue K, arithFrobenius_mem_unramLevelGeneratorSet K,
    arithFrobenius_eq_self_iff K,
    fun _π hπmax f hf0 hf1 hf _u hu g hg0 hg1 hg =>
      exists_dworkTheta_step2_of_residue K (arithFrobenius_residue K) hq hπmax f hf0 hf1 hf hu
        g hg0 hg1 hg⟩

end ABC3.Found.PGC
