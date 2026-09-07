import ABC3.Found.PGC.PrincipalUnitsRank
import ABC3.Found.PGC.AbsoluteRamification
import ABC3.Found.PGC.RamificationImageStage

/-!
# `p^{-r}·U^v_K` は p 進対数で `𝒪_K ⊆ K` に対応する —— [pGC] Proposition 2.2 の前半

原典 [pGC] 物理 p.4 末尾〜p.5 冒頭(Proposition 2.2 の証明の第一段)。
逐語(投影の 3 変種で `U v_K` の並びが揺れるので、`>` ではなくインデント引用にしてある):

    Suppose
    that we are given (in addition to ΓK) the subgroup Γv
    K ⊆ΓK for some v = r · eK, where
    r ≥2 is an integer. Then it follows that we know the subgroup U v
    K ⊆UK. Moreover,
    it follows from the theory of the p-adic logarithm (see, e.g., [5], Chapter IV, §1) that the
    submodule

(次頁に続く)

    p−r · U v
    K ⊆UK ⊗Zp Qp
    corresponds to the submodule OK ⊆K under the isomorphism induced by the p-adic
    logarithm.

★本ファイルが埋めるのは **Moreover 以降**、すなわち

```
v = r·e_K,  r ≥ 2   ⟹   p^{-r} · log(U^v_K)  =  𝒪_K      (K の中の等式)
```

である。前半(`Γ^v_K` から `U^v_K` が分かる)は
`Found/PGC/RamificationImageStage.lean::map_ramificationFiltration_reciprocityUnits_eq_principalUnits`
で既に閉じているので、§6 で両者を繋いだ形も置く。

## 数学の中身(なぜ正しいか)

`e_K` を絶対分岐指数とすると `𝔪_K^{e_K} = p·𝒪_K` である(これが `e_K` の定義そのもの)。
したがって `v = r·e_K` に対し

```
log(U^v_K) = 𝔪_K^{r·e_K} = p^r · 𝒪_K,   ゆえに   p^{-r} · log(U^v_K) = 𝒪_K.
```

`p^{-r}` を掛ける操作は、原典では `U_K ⊗_{ℤ_p} ℚ_p` の中で行われるが、
その `U_K ⊗ ℚ_p` は p 進対数で `K` と同一視される。★**本ファイルは
テンソル積を構成せず、対数を取った後の `K` の側で `p^{-r}` 倍する**
(原典自身「under the isomorphism induced by the p-adic logarithm」と言っている)。
これは逸脱の記録 1 である。

## ★★★`r ≥ 2` はどこで効くか(名指し)

2 か所で、どちらも `norm_natCast_p_pow_le_quarter`(`‖p‖^r ≤ (1/2)^r ≤ 1/4`)を経由する。

1. ★`padicLog_bijOn_ball` を半径 `‖p‖^r` の球に当てるところ。木にある
   `Found/PGC/PadicLogSurjective.lean::padicLog_bijOn` は**半径 1/4 の球でしか
   全単射ではない**ので、`𝔪_K^{r·e_K} = p^r·𝒪_K` がその球に収まっている必要がある。
   `‖p‖ ≤ 1/2` だから `r ≥ 2` で `‖p‖^r ≤ 1/4`。
2. ★`mem_integers_scaled_padicLog`(`p^{-r}·log(u-1)` が `𝒪_K` に入ること)。
   ここでも `‖log(u-1)‖ = ‖u-1‖ ≤ ‖p‖^r` を使うために `‖p‖^r ≤ 1/4` が要る。

加えて `r ≠ 0` が `r·e_K ≠ 0` に効く(`image_sub_one_principalUnits` の `hn`)。
`r = 0` だと `U^0_K = 𝒪_K^×` で対数は全単射にならない。

★**古典的な最良の条件は `v > e_K/(p-1)`** であり、`r = 1`(すなわち `v = e_K`)でも
`p ≥ 3` なら正しい。★本ファイルの証明が実際に使っているのは `‖p‖^r ≤ 1/4`
(すなわち `p^r ≥ 4`)であって、`r = 1` では `p ≥ 5` のときだけ通る。
★**原典の `r ≥ 2` は `p` に依らず一様に通る十分条件**であり、本ファイルはそれを採る。

## 構成

* §1 ★**抽象核** —— 分岐・付値・Galois・p 進の語彙が 1 語も出てこない 2 本。
  `smul_setOf_norm_le`(ノルム付き可除環での球のスカラー倍)と
  `bijOn_setOf_norm_le_of_norm_eq`(★`Norm E` だけ。ノルムを保つ球の全単射は
  小さい球でも全単射)。
* §2 絶対分岐指数 `e_K` の基本性質。★木にはこれまで
  `absoluteRamificationIndex` についての補題が `e·f = [K:ℚ_p]` 以外に
  **1 本も無かった**ので、ここで `span{p} = 𝔪^{e_K}`・`‖p‖ = ‖π‖^{e_K}`・
  `e_K ≠ 0` を入れる。
* §3 p 進対数を半径 `s ≤ 1/4` の球に制限する。
* §4 `U^n_K` の `u ↦ u-1` による像は閉球 `‖x‖ ≤ ‖π‖^n`。
* §5 主定理(集合の形と群同型の形)。
* §6 `Γ^v_K` から出発した形(Lubin-Tate のデータ込み)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★**`U_K ⊗_{ℤ_p} ℚ_p` を構成していない。** 上記のとおり、対数の像側(`K` の中)で
   `p^{-r}` 倍する形に読み替えた。原典の主張と同値だが、テンソル積の圏論的な形は無い。
2. ★**`v` は自然数**(`v = r·e_K`、`r : ℕ`)。原典の `Γ^v_K` は実数添字の上付き
   分岐群だが、`RamificationImageStage.lean` が自然数の `v` でしか像を計算していない
   ので、それに合わせた。同ファイルの「逸脱の記録 1」の引き継ぎ。
3. ★**`U^v_K` は `𝒪_K^×` の部分群 `principalUnits K π v`** として実現している。
   `Γ^ab_K` の側ではなく `𝒪_K^×` 成分で述べているのも同ファイルの引き継ぎ。
4. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ)。

## ★原典より短い道(名指し)

* ★★**`𝔪^{r·e_K} = p^r·𝒪_K` を経由せずに済む。** 原典(および教科書)は
  イデアルの言葉で `𝔪^v = p^r 𝒪_K` を出してから対数を取るが、本ファイルは
  **ノルムだけ**で通す: `‖p‖ = ‖π‖^{e_K}` さえ出れば
  `{x | ‖x‖ ≤ ‖π‖^{r e_K}} = {x | ‖x‖ ≤ ‖p‖^r}` は `pow_mul` 1 回である。
  イデアル版(`maximalIdeal_pow_mul_absoluteRamificationIndex`)は**系として**置いてある。
* ★**`e_K = ramificationIdx` の同定に Dedekind 環の因子分解が要らない。**
  `Ideal.ramificationIdx_spec` の第 2 条件(`¬ span{p} ≤ 𝔪^{e+1}`)は
  木にある `norm_le_of_mem_span_pow` と `norm_pi_lt_one` で**ノルムの不等式 1 本**に落ちる。
* ★**対数の全単射性を小さい球へ落とすのに解析は一切要らない。**
  `norm_padicLog_eq`(ノルム保存)があれば §1 の `Norm` だけの補題で機械的に落ちる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Pointwise

variable {p : ℕ} [Fact p.Prime]

/-! ## §1 抽象核

★この 2 本には分岐・付値・Galois・p 進の語彙が 1 語も出てこない。 -/

/-- ★★**抽象核 1(ノルム付き可除環)** —— `c ≠ 0` なら
`c • {x | ‖x‖ ≤ s} = {y | ‖y‖ ≤ ‖c‖·s}`。

「`p^{-r}` を掛けると閉球の半径が `‖p‖^{-r}` 倍になる」の中身はこれだけ。 -/
theorem smul_setOf_norm_le {F : Type*} [NormedDivisionRing F] {c : F} (hc : c ≠ 0) (s : ℝ) :
    c • {x : F | ‖x‖ ≤ s} = {y : F | ‖y‖ ≤ ‖c‖ * s} := by
  have hcpos : (0 : ℝ) < ‖c‖ := norm_pos_iff.mpr hc
  ext y
  constructor
  · rintro ⟨x, hx, rfl⟩
    simp only [Set.mem_setOf_eq, smul_eq_mul, norm_mul]
    exact mul_le_mul_of_nonneg_left hx (le_of_lt hcpos)
  · intro hy
    refine ⟨c⁻¹ * y, ?_, ?_⟩
    · simp only [Set.mem_setOf_eq, norm_mul, norm_inv]
      rw [inv_mul_le_iff₀ hcpos]
      exact hy
    · simp only [smul_eq_mul]
      exact mul_inv_cancel_left₀ hc y

/-- ★★**抽象核 2(`Norm E` だけ)** —— ノルムを保つ球の全単射は、
**より小さい球でも全単射**。

これが「p 進対数は半径 1/4 の球でしか全単射でない」という制約を
半径 `‖p‖^r` の球へ運ぶ唯一の道具である。代数構造も位相も要らない。 -/
theorem bijOn_setOf_norm_le_of_norm_eq {E : Type*} [Norm E] {f : E → E} {t s : ℝ}
    (hbij : Set.BijOn f {x : E | ‖x‖ ≤ t} {y : E | ‖y‖ ≤ t})
    (hnorm : ∀ x : E, ‖x‖ ≤ t → ‖f x‖ = ‖x‖) (hst : s ≤ t) :
    Set.BijOn f {x : E | ‖x‖ ≤ s} {y : E | ‖y‖ ≤ s} := by
  refine ⟨fun x hx => ?_, hbij.injOn.mono (fun x hx => le_trans hx hst), fun y hy => ?_⟩
  · simp only [Set.mem_setOf_eq] at hx ⊢
    rw [hnorm x (le_trans hx hst)]; exact hx
  · obtain ⟨x, hx, hfx⟩ := hbij.surjOn (show y ∈ {y : E | ‖y‖ ≤ t} from le_trans hy hst)
    refine ⟨x, ?_, hfx⟩
    simp only [Set.mem_setOf_eq] at hx ⊢
    rw [← hnorm x hx, hfx]; exact hy

/-! ## §2 絶対分岐指数 `e_K`

★`Found/PGC/AbsoluteRamification.lean` は `absoluteRamificationIndex` を定義したが、
`e·f = [K:ℚ_p]` 以外に補題を 1 本も持っていなかった。ここで
「`p` は素元の `e_K` 乗と同伴」という**定義そのもの**を取り出す。 -/

/-- `p ≠ 0`(`𝒪_K` の中で)。 -/
theorem natCast_p_ne_zero_integers (K : PAdicLocalField p) :
    ((p : ℕ) : 𝒪[K.carrier]) ≠ 0 := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  intro h
  have h1 : ((((p : ℕ) : 𝒪[K.carrier])) : K.carrier) = ((0 : 𝒪[K.carrier]) : K.carrier) :=
    congrArg (fun z : 𝒪[K.carrier] => (z : K.carrier)) h
  push_cast at h1
  exact (Nat.cast_ne_zero (R := K.carrier)).mpr (Nat.Prime.ne_zero Fact.out) h1

/-- `0 < ‖p‖`。 -/
theorem norm_natCast_p_pos (K : PAdicLocalField p) : (0 : ℝ) < ‖((p : ℕ) : K.carrier)‖ := by
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  exact norm_pos_iff.mpr ((Nat.cast_ne_zero (R := K.carrier)).mpr (Nat.Prime.ne_zero Fact.out))

/-- `p` は素元 `π` の何乗かと同伴(`𝒪_K` は DVR)。 -/
theorem exists_uniformizer_pow_eq_p (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    ∃ e : ℕ, ‖((p : ℕ) : K.carrier)‖ = ‖(π : K.carrier)‖ ^ e ∧
      Ideal.span ({((p : ℕ) : 𝒪[K.carrier])} : Set 𝒪[K.carrier]) = (Ideal.span {π}) ^ e := by
  haveI := valuationRing_isDVR K
  have hirr : Irreducible π := (IsDiscreteValuationRing.irreducible_iff_uniformizer π).2 hπmax
  obtain ⟨e, u, hu⟩ :=
    IsDiscreteValuationRing.eq_unit_mul_pow_irreducible (natCast_p_ne_zero_integers K) hirr
  refine ⟨e, ?_, ?_⟩
  · have h1 : (((p : ℕ) : 𝒪[K.carrier]) : K.carrier)
        = ((u : 𝒪[K.carrier]) : K.carrier) * ((π : K.carrier)) ^ e := by
      rw [hu]; push_cast; ring
    have hun : ‖((u : 𝒪[K.carrier]) : K.carrier)‖ = 1 :=
      Valued.integer.isUnit_iff_norm_eq_one.mp u.isUnit
    have h3 : (((p : ℕ) : 𝒪[K.carrier]) : K.carrier) = ((p : ℕ) : K.carrier) := by push_cast; ring
    rw [← h3, h1, norm_mul, hun, one_mul, norm_pow]
  · rw [Ideal.span_singleton_pow]
    exact (Ideal.span_singleton_eq_span_singleton.mpr
      (⟨u, by rw [hu]; ring⟩ : Associated (π ^ e) ((p : ℕ) : 𝒪[K.carrier]))).symm

/-- ★**`e_K` の同定** —— `span{p} = 𝔪^e` なら `absoluteRamificationIndex K = e`。

`Ideal.ramificationIdx_spec` の第 2 条件は、木にある `norm_le_of_mem_span_pow` と
`norm_pi_lt_one` で**ノルムの不等式 1 本**に落ちる(Dedekind 環の因子分解は不要)。 -/
theorem absoluteRamificationIndex_eq_of_span (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) {e : ℕ}
    (hnorm : ‖((p : ℕ) : K.carrier)‖ = ‖(π : K.carrier)‖ ^ e)
    (he : Ideal.span ({((p : ℕ) : 𝒪[K.carrier])} : Set 𝒪[K.carrier]) = (Ideal.span {π}) ^ e) :
    absoluteRamificationIndex K = e := by
  haveI := valuationRing_isDVR K
  have hirr : Irreducible π := (IsDiscreteValuationRing.irreducible_iff_uniformizer π).2 hπmax
  have hπ0 : (π : K.carrier) ≠ 0 := fun h => hirr.ne_zero (Subtype.ext h)
  have hπpos : (0 : ℝ) < ‖(π : K.carrier)‖ := norm_pos_iff.mpr hπ0
  have hmap : Ideal.map (algebraMap ℤ_[p] 𝒪[K.carrier]) (IsLocalRing.maximalIdeal ℤ_[p])
      = Ideal.span ({((p : ℕ) : 𝒪[K.carrier])} : Set 𝒪[K.carrier]) := by
    rw [PadicInt.maximalIdeal_eq_span_p, Ideal.map_span]; simp
  apply Ideal.ramificationIdx_spec
  · rw [hmap, hπmax, he]
  · rw [hmap, hπmax]
    intro hcon
    have hmem : ((p : ℕ) : 𝒪[K.carrier]) ∈ (Ideal.span ({π} : Set 𝒪[K.carrier])) ^ (e + 1) :=
      hcon (Ideal.mem_span_singleton_self _)
    rw [Ideal.span_singleton_pow] at hmem
    have hn := norm_le_of_mem_span_pow K (e + 1) _ hmem
    have h3 : (((p : ℕ) : 𝒪[K.carrier]) : K.carrier) = ((p : ℕ) : K.carrier) := by push_cast; ring
    rw [h3, hnorm, pow_succ] at hn
    nlinarith [pow_pos hπpos e, norm_pi_lt_one K hπmax]

/-- ★**`‖p‖ = ‖π‖^{e_K}`** —— 絶対分岐指数のノルムによる特徴づけ。 -/
theorem norm_natCast_p_eq_norm_uniformizer_pow (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) :
    ‖((p : ℕ) : K.carrier)‖ = ‖(π : K.carrier)‖ ^ (absoluteRamificationIndex K) := by
  obtain ⟨e, hnorm, hspan⟩ := exists_uniformizer_pow_eq_p K hπmax
  rw [absoluteRamificationIndex_eq_of_span K hπmax hnorm hspan]
  exact hnorm

/-- ★**`p·𝒪_K = 𝔪_K^{e_K}`** —— 絶対分岐指数のイデアルによる特徴づけ
(素元の取り方に依らない形)。 -/
theorem span_natCast_p_eq_maximalIdeal_pow (K : PAdicLocalField p) :
    Ideal.span ({((p : ℕ) : 𝒪[K.carrier])} : Set 𝒪[K.carrier])
      = (IsLocalRing.maximalIdeal 𝒪[K.carrier]) ^ (absoluteRamificationIndex K) := by
  haveI := valuationRing_isDVR K
  obtain ⟨π, hirr⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  have hπmax := (IsDiscreteValuationRing.irreducible_iff_uniformizer π).1 hirr
  obtain ⟨e, hnorm, hspan⟩ := exists_uniformizer_pow_eq_p K hπmax
  rw [hπmax, absoluteRamificationIndex_eq_of_span K hπmax hnorm hspan]
  exact hspan

/-- `e_K ≠ 0`(`e_K = 0` なら `‖p‖ = 1` になり `‖p‖ ≤ 1/2` に反する)。 -/
theorem absoluteRamificationIndex_ne_zero (K : PAdicLocalField p) :
    absoluteRamificationIndex K ≠ 0 := by
  haveI := valuationRing_isDVR K
  obtain ⟨π, hirr⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  have hπmax := (IsDiscreteValuationRing.irreducible_iff_uniformizer π).1 hirr
  intro h
  have hn := norm_natCast_p_eq_norm_uniformizer_pow K hπmax
  rw [h, pow_zero] at hn
  have h2 := norm_natCast_p_le_half K
  rw [hn] at h2
  norm_num at h2

/-- ★**`𝔪_K^{r·e_K} = p^r·𝒪_K`** —— イデアル版(原典が経由する形)。
★本ファイルの主定理はこれを**使わない**(ノルムだけで通す)。系として置く。 -/
theorem maximalIdeal_pow_mul_absoluteRamificationIndex (K : PAdicLocalField p) (r : ℕ) :
    (IsLocalRing.maximalIdeal 𝒪[K.carrier]) ^ (r * absoluteRamificationIndex K)
      = Ideal.span ({((p : ℕ) : 𝒪[K.carrier]) ^ r} : Set 𝒪[K.carrier]) := by
  rw [pow_mul', ← span_natCast_p_eq_maximalIdeal_pow, Ideal.span_singleton_pow]

/-- ★**`‖π‖^{r·e_K} = ‖p‖^r`** —— 本ファイルが実際に使う形。 -/
theorem norm_uniformizer_pow_mul_ramificationIndex (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) (r : ℕ) :
    ‖(π : K.carrier)‖ ^ (r * absoluteRamificationIndex K) = ‖((p : ℕ) : K.carrier)‖ ^ r := by
  rw [norm_natCast_p_eq_norm_uniformizer_pow K hπmax, ← pow_mul, mul_comm]

/-- ★★**`r ≥ 2` が効く唯一の場所** —— `‖p‖ ≤ 1/2` だから `‖p‖^r ≤ (1/2)^2 = 1/4`。 -/
theorem norm_natCast_p_pow_le_quarter (K : PAdicLocalField p) {r : ℕ} (hr : 2 ≤ r) :
    ‖((p : ℕ) : K.carrier)‖ ^ r ≤ 1 / 4 := by
  calc ‖((p : ℕ) : K.carrier)‖ ^ r ≤ (1 / 2 : ℝ) ^ r :=
        pow_le_pow_left₀ (norm_nonneg _) (norm_natCast_p_le_half K) r
    _ ≤ (1 / 2 : ℝ) ^ 2 := pow_le_pow_of_le_one (by norm_num) (by norm_num) hr
    _ = 1 / 4 := by norm_num

/-! ## §3 p 進対数を小さい球へ制限する -/

/-- `log 0 = 0`(`padicLog_mul` に `x = y = 0` を入れる)。 -/
theorem padicLog_zero (K : PAdicLocalField p) : padicLog K 0 = 0 := by
  have h2 := padicLog_mul (K := K) (x := 0) (y := 0) (by simp) (by simp)
  simpa using h2

/-- `Found/PGC/PadicLogInjective.lean::norm_padicLog_eq` の `x = 0` を込めた形。 -/
theorem norm_padicLog_eq' (K : PAdicLocalField p) {x : K.carrier} (hx : ‖x‖ ≤ 1 / 4) :
    ‖padicLog K x‖ = ‖x‖ := by
  rcases eq_or_ne x 0 with rfl | h
  · simp [padicLog_zero]
  · exact norm_padicLog_eq K hx h

/-- ★**p 進対数は半径 `s ≤ 1/4` のどの閉球でも全単射**。
`padicLog_bijOn`(半径 1/4)＋ノルム保存＋§1 の抽象核。 -/
theorem padicLog_bijOn_ball (K : PAdicLocalField p) {s : ℝ} (hs : s ≤ 1 / 4) :
    Set.BijOn (padicLog K) {x : K.carrier | ‖x‖ ≤ s} {y : K.carrier | ‖y‖ ≤ s} :=
  bijOn_setOf_norm_le_of_norm_eq (padicLog_bijOn K) (fun _ hx => norm_padicLog_eq' K hx) hs

/-! ## §4 `U^n_K` の像 -/

/-- `u ∈ U^n_K` なら `‖u-1‖ ≤ ‖π‖^n`。 -/
theorem norm_sub_one_le_of_mem_principalUnits (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    {n : ℕ} {u : (𝒪[K.carrier])ˣ} (hu : u ∈ principalUnits K π n) :
    ‖((u : 𝒪[K.carrier]) : K.carrier) - 1‖ ≤ ‖(π : K.carrier)‖ ^ n := by
  have hmem : (u : 𝒪[K.carrier]) - 1 ∈ Ideal.span ({π ^ n} : Set (𝒪[K.carrier])) := hu
  have h := norm_le_of_mem_span_pow K n _ hmem
  simpa using h

/-- ★**`U^n_K - 1 = 𝔪_K^n`**(`n ≠ 0`)—— `u ↦ u-1` による `U^n_K` の像は
ちょうど閉球 `‖x‖ ≤ ‖π‖^n`。逆向きは `‖x‖ < 1` から `‖1+x‖ = 1`、
すなわち `1+x` が `𝒪_K` の単数であることによる。 -/
theorem image_sub_one_principalUnits (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) (hπne0 : π ≠ 0)
    {n : ℕ} (hn : n ≠ 0) :
    (fun u : (𝒪[K.carrier])ˣ => ((u : 𝒪[K.carrier]) : K.carrier) - 1) ''
        (principalUnits K π n : Set (𝒪[K.carrier])ˣ)
      = {x : K.carrier | ‖x‖ ≤ ‖(π : K.carrier)‖ ^ n} := by
  have hlt : ‖(π : K.carrier)‖ ^ n < 1 :=
    pow_lt_one₀ (norm_nonneg _) (norm_pi_lt_one K hπmax) hn
  ext x
  constructor
  · rintro ⟨u, hu, rfl⟩
    exact norm_sub_one_le_of_mem_principalUnits K hu
  · intro hx
    simp only [Set.mem_setOf_eq] at hx
    have hxlt : ‖x‖ < 1 := lt_of_le_of_lt hx hlt
    have hone : ‖(1 : K.carrier) + x‖ = 1 := by
      rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier)]
      · rw [norm_one]; exact max_eq_left (le_of_lt hxlt)
      · rw [norm_one]; exact ne_of_gt hxlt
    have hmemO : (1 : K.carrier) + x ∈ 𝒪[K.carrier] := by
      rw [Valued.integer.mem_iff, hone]
    set y : 𝒪[K.carrier] := ⟨(1 : K.carrier) + x, hmemO⟩ with hy
    have hyunit : IsUnit y := Valued.integer.isUnit_iff_norm_eq_one.mpr hone
    refine ⟨hyunit.unit, ?_, ?_⟩
    · show ((hyunit.unit : 𝒪[K.carrier]) : 𝒪[K.carrier]) - 1 ∈
        Ideal.span ({π ^ n} : Set (𝒪[K.carrier]))
      rw [hyunit.unit_spec]
      apply mem_span_pow_of_norm_le K hπne0
      have hval : ((y - 1 : 𝒪[K.carrier]) : K.carrier) = x := by simp [hy]
      rw [hval]; exact hx
    · show ((hyunit.unit : 𝒪[K.carrier]) : K.carrier) - 1 = x
      rw [hyunit.unit_spec]
      simp [hy]

/-! ## §5 ★★★主定理 -/

/-- ★★★★**[pGC] Proposition 2.2 の第一段(集合の形)** ——

    p^{-r} · log(U^v_K)  =  𝒪_K        (v = r·e_K,  r ≥ 2)

原典の「the submodule `p^{-r}·U^v_K ⊆ U_K ⊗_{ℤ_p} ℚ_p` corresponds to the submodule
`𝒪_K ⊆ K` under the isomorphism induced by the p-adic logarithm」。
★テンソル積は構成せず、対数を取った後の `K` の中で `p^{-r}` 倍している(逸脱 1)。 -/
theorem smul_padicLog_image_principalUnits_eq_integers (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0) {r : ℕ} (hr : 2 ≤ r) :
    (((p : ℕ) : K.carrier) ^ r)⁻¹ •
        (padicLog K '' ((fun u : (𝒪[K.carrier])ˣ => ((u : 𝒪[K.carrier]) : K.carrier) - 1) ''
          (principalUnits K π (r * absoluteRamificationIndex K) : Set (𝒪[K.carrier])ˣ)))
      = (𝒪[K.carrier] : Set K.carrier) := by
  have hne : r * absoluteRamificationIndex K ≠ 0 :=
    Nat.mul_ne_zero (by omega) (absoluteRamificationIndex_ne_zero K)
  have hppos : (0 : ℝ) < ‖((p : ℕ) : K.carrier)‖ ^ r := pow_pos (norm_natCast_p_pos K) r
  have hpne : (((p : ℕ) : K.carrier) ^ r)⁻¹ ≠ 0 := by
    refine inv_ne_zero (pow_ne_zero r ?_)
    exact norm_pos_iff.mp (norm_natCast_p_pos K)
  rw [image_sub_one_principalUnits K hπmax hπne0 hne,
    norm_uniformizer_pow_mul_ramificationIndex K hπmax r,
    (padicLog_bijOn_ball K (norm_natCast_p_pow_le_quarter K hr)).image_eq,
    smul_setOf_norm_le hpne]
  have hone : ‖(((p : ℕ) : K.carrier) ^ r)⁻¹‖ * ‖((p : ℕ) : K.carrier)‖ ^ r = 1 := by
    rw [norm_inv, norm_pow]
    exact inv_mul_cancel₀ (ne_of_gt hppos)
  rw [hone]
  ext y
  exact (Valued.integer.mem_iff).symm

/-- 原典が経由する中間形 `log(U^v_K) = p^r·𝒪_K`。 -/
theorem padicLog_image_principalUnits_eq_smul_integers (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0) {r : ℕ} (hr : 2 ≤ r) :
    padicLog K '' ((fun u : (𝒪[K.carrier])ˣ => ((u : 𝒪[K.carrier]) : K.carrier) - 1) ''
        (principalUnits K π (r * absoluteRamificationIndex K) : Set (𝒪[K.carrier])ˣ))
      = ((p : ℕ) : K.carrier) ^ r • (𝒪[K.carrier] : Set K.carrier) := by
  have hne : r * absoluteRamificationIndex K ≠ 0 :=
    Nat.mul_ne_zero (by omega) (absoluteRamificationIndex_ne_zero K)
  have hpne : ((p : ℕ) : K.carrier) ^ r ≠ 0 :=
    pow_ne_zero r (norm_pos_iff.mp (norm_natCast_p_pos K))
  have hOK : (𝒪[K.carrier] : Set K.carrier) = {y : K.carrier | ‖y‖ ≤ 1} := by
    ext y; exact Valued.integer.mem_iff
  rw [image_sub_one_principalUnits K hπmax hπne0 hne,
    norm_uniformizer_pow_mul_ramificationIndex K hπmax r,
    (padicLog_bijOn_ball K (norm_natCast_p_pow_le_quarter K hr)).image_eq,
    hOK, smul_setOf_norm_le hpne, norm_pow, mul_one]

/-- `p^{-r}·log(u-1)` は `𝒪_K` に入る(`r ≥ 2` が効く 2 か所目)。 -/
theorem mem_integers_scaled_padicLog (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) {r : ℕ} (hr : 2 ≤ r)
    {u : (𝒪[K.carrier])ˣ} (hu : u ∈ principalUnits K π (r * absoluteRamificationIndex K)) :
    (((p : ℕ) : K.carrier) ^ r)⁻¹ * padicLog K (((u : 𝒪[K.carrier]) : K.carrier) - 1)
      ∈ 𝒪[K.carrier] := by
  have hppos : (0 : ℝ) < ‖((p : ℕ) : K.carrier)‖ ^ r := pow_pos (norm_natCast_p_pos K) r
  have hb : ‖((u : 𝒪[K.carrier]) : K.carrier) - 1‖ ≤ ‖((p : ℕ) : K.carrier)‖ ^ r := by
    have h := norm_sub_one_le_of_mem_principalUnits K hu
    rwa [norm_uniformizer_pow_mul_ramificationIndex K hπmax r] at h
  have hlog : ‖padicLog K (((u : 𝒪[K.carrier]) : K.carrier) - 1)‖
      ≤ ‖((p : ℕ) : K.carrier)‖ ^ r := by
    rw [norm_padicLog_eq' K (le_trans hb (norm_natCast_p_pow_le_quarter K hr))]
    exact hb
  rw [Valued.integer.mem_iff, norm_mul, norm_inv, norm_pow]
  rw [inv_mul_le_iff₀ hppos, mul_one]
  exact hlog

/-- ★**`u ↦ p^{-r}·log(u-1)` は群準同型** `U^v_K →* Multiplicative 𝒪_K`。
乗法から加法への変換は `padicLog_mul`(`log((1+x)(1+y)) = log(1+x)+log(1+y)`)。 -/
noncomputable def padicLogPrincipalUnitsHom (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) {r : ℕ} (hr : 2 ≤ r) :
    principalUnits K π (r * absoluteRamificationIndex K) →* Multiplicative 𝒪[K.carrier] where
  toFun u := Multiplicative.ofAdd
    (⟨(((p : ℕ) : K.carrier) ^ r)⁻¹ *
        padicLog K ((((u : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1),
      mem_integers_scaled_padicLog K hπmax hr u.2⟩ : 𝒪[K.carrier])
  map_one' := by
    apply congrArg Multiplicative.ofAdd
    apply Subtype.ext
    show (((p : ℕ) : K.carrier) ^ r)⁻¹ * padicLog K
      (((((1 : principalUnits K π (r * absoluteRamificationIndex K)) :
        (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1) = 0
    have h : (((((1 : principalUnits K π (r * absoluteRamificationIndex K)) :
        (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1) = 0 := by simp
    rw [h, padicLog_zero, mul_zero]
  map_mul' a b := by
    apply congrArg Multiplicative.ofAdd
    apply Subtype.ext
    show (((p : ℕ) : K.carrier) ^ r)⁻¹ * padicLog K
        ((((((a * b : principalUnits K π (r * absoluteRamificationIndex K)) :
          (𝒪[K.carrier])ˣ)) : 𝒪[K.carrier]) : K.carrier) - 1)
      = (((p : ℕ) : K.carrier) ^ r)⁻¹ * padicLog K
          ((((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
        + (((p : ℕ) : K.carrier) ^ r)⁻¹ * padicLog K
          ((((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
    have hq : ‖((p : ℕ) : K.carrier)‖ ^ r ≤ 1 / 4 := norm_natCast_p_pow_le_quarter K hr
    have ha : ‖(((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1‖ < 1 := by
      have h := norm_sub_one_le_of_mem_principalUnits K a.2
      rw [norm_uniformizer_pow_mul_ramificationIndex K hπmax r] at h
      linarith [le_trans h hq]
    have hb : ‖(((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1‖ < 1 := by
      have h := norm_sub_one_le_of_mem_principalUnits K b.2
      rw [norm_uniformizer_pow_mul_ramificationIndex K hπmax r] at h
      linarith [le_trans h hq]
    have hexp : (((((a * b : principalUnits K π (r * absoluteRamificationIndex K)) :
          (𝒪[K.carrier])ˣ)) : 𝒪[K.carrier]) : K.carrier) - 1
        = ((((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
          + ((((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
          + ((((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
            * ((((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1) := by
      push_cast; ring
    rw [hexp, padicLog_mul ha hb, mul_add]

/-- 単射は `padicLog_injOn`、全射は `padicLog_bijOn_ball` ＋ `image_sub_one_principalUnits`。 -/
theorem padicLogPrincipalUnitsHom_bijective (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) (hπne0 : π ≠ 0)
    {r : ℕ} (hr : 2 ≤ r) :
    Function.Bijective (padicLogPrincipalUnitsHom K hπmax hr) := by
  have hq : ‖((p : ℕ) : K.carrier)‖ ^ r ≤ 1 / 4 := norm_natCast_p_pow_le_quarter K hr
  have hpne : ((p : ℕ) : K.carrier) ^ r ≠ 0 :=
    pow_ne_zero r (norm_pos_iff.mp (norm_natCast_p_pos K))
  have hne : r * absoluteRamificationIndex K ≠ 0 :=
    Nat.mul_ne_zero (by omega) (absoluteRamificationIndex_ne_zero K)
  constructor
  · intro a b hab
    have h0 : (((p : ℕ) : K.carrier) ^ r)⁻¹ *
          padicLog K ((((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
        = (((p : ℕ) : K.carrier) ^ r)⁻¹ *
          padicLog K ((((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1) :=
      congrArg Subtype.val hab
    have h1 : padicLog K ((((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1)
        = padicLog K ((((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1) :=
      mul_left_cancel₀ (inv_ne_zero hpne) h0
    have ha : (((a : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1
        ∈ {x : K.carrier | ‖x‖ ≤ 1 / 4} := by
      have h := norm_sub_one_le_of_mem_principalUnits K a.2
      rw [norm_uniformizer_pow_mul_ramificationIndex K hπmax r] at h
      exact le_trans h hq
    have hb : (((b : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1
        ∈ {x : K.carrier | ‖x‖ ≤ 1 / 4} := by
      have h := norm_sub_one_le_of_mem_principalUnits K b.2
      rw [norm_uniformizer_pow_mul_ramificationIndex K hπmax r] at h
      exact le_trans h hq
    have hsub := padicLog_injOn K ha hb h1
    apply Subtype.ext
    apply Units.ext
    apply Subtype.ext
    linear_combination hsub
  · intro y
    set w : K.carrier :=
      ((p : ℕ) : K.carrier) ^ r * ((Multiplicative.toAdd y : 𝒪[K.carrier]) : K.carrier) with hw
    have hyle : ‖((Multiplicative.toAdd y : 𝒪[K.carrier]) : K.carrier)‖ ≤ 1 :=
      Valued.integer.mem_iff.mp (Multiplicative.toAdd y).2
    have hwle : ‖w‖ ≤ ‖((p : ℕ) : K.carrier)‖ ^ r := by
      rw [hw, norm_mul, norm_pow]
      nlinarith [pow_nonneg (norm_nonneg (((p : ℕ) : K.carrier))) r,
        norm_nonneg (((Multiplicative.toAdd y : 𝒪[K.carrier]) : K.carrier))]
    obtain ⟨x, hx, hlogx⟩ := (padicLog_bijOn_ball K hq).surjOn hwle
    have hxmem : x ∈ {x : K.carrier |
        ‖x‖ ≤ ‖(π : K.carrier)‖ ^ (r * absoluteRamificationIndex K)} := by
      rw [norm_uniformizer_pow_mul_ramificationIndex K hπmax r]; exact hx
    rw [← image_sub_one_principalUnits K hπmax hπne0 hne] at hxmem
    obtain ⟨u, hu, hux⟩ := hxmem
    have hux' : (((u : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1 = x := hux
    refine ⟨⟨u, hu⟩, ?_⟩
    apply congrArg Multiplicative.ofAdd
    apply Subtype.ext
    show (((p : ℕ) : K.carrier) ^ r)⁻¹ * padicLog K (((u : 𝒪[K.carrier]) : K.carrier) - 1)
      = ((Multiplicative.toAdd y : 𝒪[K.carrier]) : K.carrier)
    rw [hux', hlogx, hw, ← mul_assoc, inv_mul_cancel₀ hpne, one_mul]

/-- ★★★★**[pGC] Proposition 2.2 の第一段(同型の形)** ——

    U^v_K  ≃*  Multiplicative 𝒪_K        (v = r·e_K,  r ≥ 2)

同型を与えているのは `u ↦ p^{-r}·log(u-1)`。★これが原典の
「corresponds to the submodule `𝒪_K ⊆ K` under the isomorphism induced by
the p-adic logarithm」の、群同型としての形である。 -/
noncomputable def padicLogPrincipalUnitsEquiv (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) (hπne0 : π ≠ 0)
    {r : ℕ} (hr : 2 ≤ r) :
    principalUnits K π (r * absoluteRamificationIndex K) ≃* Multiplicative 𝒪[K.carrier] :=
  MulEquiv.ofBijective (padicLogPrincipalUnitsHom K hπmax hr)
    (padicLogPrincipalUnitsHom_bijective K hπmax hπne0 hr)

def padicLogPrincipalUnitsEquiv.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- 同型の値の式(`rfl`)。 -/
theorem padicLogPrincipalUnitsEquiv_apply (K : PAdicLocalField p) {π : 𝒪[K.carrier]}
    (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π}) (hπne0 : π ≠ 0)
    {r : ℕ} (hr : 2 ≤ r) (u : principalUnits K π (r * absoluteRamificationIndex K)) :
    ((Multiplicative.toAdd (padicLogPrincipalUnitsEquiv K hπmax hπne0 hr u) :
        𝒪[K.carrier]) : K.carrier)
      = (((p : ℕ) : K.carrier) ^ r)⁻¹ *
        padicLog K ((((u : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]) : K.carrier) - 1) := rfl

/-! ## §6 ★★`Γ^v_K` から出発した形

`Found/PGC/RamificationImageStage.lean` の
`map_ramificationFiltration_reciprocityUnits_eq_principalUnits`(`Art(Γ^n_K) = U^n_K`)
と繋ぐと、原典の段落全体——「`Γ^v_K` が分かれば `U^v_K` が分かり、`p^{-r}·U^v_K` は
`𝒪_K` に対応する」——が 1 本になる。 -/

section LubinTate

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal 𝒪[K.carrier]) 𝒪[K.carrier]]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField 𝒪[K.carrier]) pp]
    [Fintype (IsLocalRing.ResidueField 𝒪[K.carrier])]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField 𝒪[K.carrier]) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal 𝒪[K.carrier] = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries 𝒪[K.carrier]) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue 𝒪[K.carrier]) f = PowerSeries.X ^ (pp ^ ff))

/-- ★★★★★**原典の段落全体** ——

    v = r·e_K,  r ≥ 2   ⟹   p^{-r} · log(Art(Γ^v_K))  =  𝒪_K.

前半 `Art(Γ^v_K) = U^v_K` は `RamificationImageStage.lean`(仮定ゼロ)、
後半 `p^{-r}·log(U^v_K) = 𝒪_K` は本ファイル §5。 -/
theorem smul_padicLog_image_ramificationFiltration_eq_integers {r : ℕ} (hr : 2 ≤ r) :
    (((p : ℕ) : K.carrier) ^ r)⁻¹ •
        (padicLog K '' ((fun u : (𝒪[K.carrier])ˣ => ((u : 𝒪[K.carrier]) : K.carrier) - 1) ''
          ((Subgroup.map (reciprocityUnits K hq hπmax hπne0 f hf0 hf1 hf)
              ((ramificationFiltration p).Gv K
                ((((r * absoluteRamificationIndex K : ℕ)) : ℕ) : ℝ)) :
            Subgroup (𝒪[K.carrier])ˣ) : Set (𝒪[K.carrier])ˣ)))
      = (𝒪[K.carrier] : Set K.carrier) := by
  rw [map_ramificationFiltration_reciprocityUnits_eq_principalUnits K hq hπmax hπne0
    f hf0 hf1 hf (r * absoluteRamificationIndex K)]
  exact smul_padicLog_image_principalUnits_eq_integers K hπmax hπne0 hr

def smul_padicLog_image_ramificationFiltration_eq_integers.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

end LubinTate

end ABC3.Found.PGC
