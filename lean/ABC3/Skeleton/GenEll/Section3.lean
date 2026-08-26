import ABC3.Meta.Claim
import ABC3.Interface.GenEll.TateLocal
import ABC3.Interface.GenEll.EllModuli
import ABC3.Found.GenEll.Lemma31
import ABC3.Found.GenEll.Elementary
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Found.GenEll.BDClass
import Mathlib.Topology.Algebra.Ring.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# [GenEll] §3 Full Special Linear Galois Actions —— 残り 8 件(`Skeleton`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.13–p.18。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.13):
> Lemma 3.1. (The Structure of SL2) Let l ≥ 5 be a prime number. Then:

`Theorem 3.8` は `Skeleton/GenEll/GaloisImage.lean` にある。ここでは残り 8 件を固定する。

## ★★このページで実測した最も重要なこと

★**`Lemma 3.6`(An Elementary Estimate)は純粋な実解析である。**
Arakelov 理論も Galois 表現も Tate 曲線も要らない——
「`ϵ > 0` に対し定数 `C₀ > 0` が存在して、`y ≥ 1` かつ `x ≥ C₀y^{1+ϵ}` なる全ての
実数 `x, y` について `x ≥ y·log(x)`」だけである。
★**したがって 2 理論なしで届く項目は 3 件ではなく 4 件**である
(`Lemma 3.1`・`Lemma 3.6`・`Lemma 4.1`・`Lemma 4.2`)。
以前の見積もりを 1 件上方修正する。

## ★`Lemma 3.1` の位置づけ

★(i)(ii)(iii) は **`Found/GenEll/Lemma31.lean` に実装済み**(`sorry` 無し)。
残るのは **(iv)** だけであり、原文はそこで **[Serre] Chapter IV, §3.4, Lemma 3** を引く。
★その [Serre] は `0_Source` に無いので、**引用ではなく自分で証明する**ことになる。

★ここで `Lemma 3.1` **全体**を 1 つの statement として固定するのは、
「(i)(ii)(iii) が済んでいる」ことと「命題全体が済んでいる」ことを
**混同しないため**である(`.src` の 2 値規則)。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Interface.GenEll ABC3.Found.GenEll
open Matrix Matrix.SpecialLinearGroup
open scoped MatrixGroups

/-! ## Lemma 3.1 —— 全体(i)–(iv) -/

/-! ★`glRedPadic`(`GL₂(ℤ_l) → GL₂(𝔽_l)` の還元)は
`Found/GenEll/Sl2Padic.lean` に移した——`Lemma 3.1, (iv)` の証明が使うからである。 -/

/-- **[GenEll] Lemma 3.1**(The Structure of SL2)。

原文 (GenEll p.13):
> Lemma 3.1. (The Structure of SL2) Let l ≥ 5 be a prime number. Then:

(i) `SL₂(𝔽_l)` は `α = (1 1 / 0 1)` と `β = (1 0 / 1 1)` で生成される。
(ii) `SL₂(𝔽_l)` の正規部分群で商が可換なものは全体。
(iii) `GL₂(𝔽_l)` の部分群が `α` と**非上三角**な行列を含めば `SL₂(𝔽_l)` を含む。
(iv) `GL₂(ℤ_l)` の**閉**部分群 `J` の `GL₂(𝔽_l)` での像が同じ条件を満たせば `SL₂(ℤ_l) ⊆ J`。

★★**4 条すべてが実装された(2026-08-17)。**
(i)(ii)(iii) は `Found/GenEll/Lemma31.lean`、
**(iv) は `Found/GenEll/Sl2Padic.lean`** ——原文が引く [Serre] Chapter IV, §3.4, Lemma 3 が
`0_Source` に無いので、段 (A)(`Sl2Level.lean`、有限群論)と
段 (B)(`Sl2Padic.lean`、位相)で**自分で証明した**。

★**本 statement はもう `sorry` ではない。**

★**「閉」を落としてはいけない**。落とせば主張が強くなる(G5 の裏返し)。
`ℤ_l` の位相から `GL₂(ℤ_l)` に入る位相で `IsClosed` を要求している。 -/
theorem lemma_3_1 (l : ℕ) [Fact (Nat.Prime l)] (hl : 5 ≤ l) :
    (Subgroup.closure ({upper (1 : ZMod l), lower (1 : ZMod l)} : Set SL(2, ZMod l)) = ⊤)
  ∧ (∀ (N : Subgroup SL(2, ZMod l)) [N.Normal],
        (∀ x y : SL(2, ZMod l) ⧸ N, x * y = y * x) → N = ⊤)
  ∧ (∀ H : Subgroup (GL (Fin 2) (ZMod l)),
        (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ H →
        (∃ M ∈ H, (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0) →
        ∀ g : SL(2, ZMod l), (toGL g : GL (Fin 2) (ZMod l)) ∈ H)
  ∧ (∀ J : Subgroup (GL (Fin 2) ℤ_[l]), IsClosed (J : Set (GL (Fin 2) ℤ_[l])) →
        (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l)) ∈ J.map (glRedPadic l) →
        (∃ M ∈ J.map (glRedPadic l), (M : Matrix (Fin 2) (Fin 2) (ZMod l)) 1 0 ≠ 0) →
        ∀ g : SL(2, ℤ_[l]), (toGL g : GL (Fin 2) ℤ_[l]) ∈ J) :=
  ABC3.Found.GenEll.lemma_3_1 l hl

/-! ## Lemma 3.2 —— l-捩れの局所階数 1 部分群 -/

/-- **[GenEll] Lemma 3.2**(Local Rank One Subgroups of l-Torsion)。

原文 (GenEll p.15):
> Lemma 3.2. (Local Rank One Subgroups of l-Torsion)

(i) `N ⊆ M_l(E)` を `G_K` 安定な 1 次元 `𝔽_l`-部分空間とすると、
**`v_K(q_E) ∈ l·ℤ` であるか、または `N = 𝔽_l(1)`**。
(ii) `E′ ≝ E/μ_l` の Tate 母数は `q_{E′} = q_E^l`、したがって
**`deg_∞(E′) = l·deg_∞(E)`**。

★(i) の「either … or」は**排他的ではない**(原文 "either … or")。 -/
theorem lemma_3_2 (D : TateLocalData) :
    (∀ (K : D.LocalField) (E : D.Curve K) (l : ℕ), Nat.Prime l →
        ∀ (N : D.StableLine K E l), l ∣ D.vq K E ∨ D.IsCyclotomic N)
  ∧ (∀ (K : D.LocalField) (E : D.Curve K) (l : ℕ), Nat.Prime l →
        D.degInf K (D.quotMu K E l) = (l : ℝ) * D.degInf K E) := by
  refine ⟨D.stableLine_dvd_or_cyclotomic, fun K E l hl => ?_⟩
  rw [D.degInf_eq K (D.quotMu K E l), D.degInf_eq K E, D.vq_quotMu K E l hl]
  push_cast
  ring

/-! ## Definition 3.3 —— 局所高さ -/

/-- **[GenEll] Definition 3.3**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★**1 文の定義**である。`Theorem 3.8` の条件 (b)「`l` is prime to the **local heights**」も、
`Corollary 4.3/4.4` の (a) も、すべてこの語を指す。 -/
def localHeight (D : TateLocalData) (K : D.LocalField) (E : D.Curve K) : ℕ :=
  D.vq K E

/-- ★原文「the **positive** integer」——正であることは定義の一部である。
`Interface` が `vq_pos` として持っているので、ここは**証明できる**。 -/
theorem localHeight_pos (D : TateLocalData) (K : D.LocalField) (E : D.Curve K) :
    0 < localHeight D K E :=
  D.vq_pos K E

/-! ## Remark 3.3.1 —— 潜在的乗法還元への拡張 -/

/-- **[GenEll] Remark 3.3.1** —— 局所高さの `ℚ` への拡張。

原文 (GenEll p.16):
> Remark 3.3.1. Note that even if EK only has potentially multiplicative reduction, one may define the local height of EK as the element ∈ Q computed by dividing the local height of

`E_K` が**潜在的**乗法還元しか持たない場合でも、乗法還元を持つ有限拡大 `L/K` を取り、
`E_K ×_K L` の局所高さを `L/K` の**分岐指数で割れば** `ℚ` の元として定義できる。 -/
noncomputable def potLocalHeight (D : TateLocalData) {K : D.LocalField} {E : D.Curve K}
    (L : D.MultExt K E) : ℚ :=
  (D.vq (D.extField L) (D.baseChange L) : ℚ) / (D.ramIdx L : ℚ)

/-- **[GenEll] Remark 3.3.1** の主張本体 —— **`L` の取り方に依らない**。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

★★**これが Remark の内容そのものである。** 上の `potLocalHeight` は
`L` を引数に取る「候補」にすぎず、`L` に依らないことを示して初めて
「`E_K` の局所高さ」と呼べる。★原文は「one verifies immediately」で済ませている。 -/
theorem potLocalHeight_indep (D : TateLocalData) {K : D.LocalField} {E : D.Curve K}
    (L L' : D.MultExt K E) :
    potLocalHeight D L = potLocalHeight D L' := by
  have he : ((D.ramIdx L : ℚ)) ≠ 0 := by
    have h := D.ramIdx_pos L
    have h0 : D.ramIdx L ≠ 0 := Nat.pos_iff_ne_zero.1 h
    exact_mod_cast h0
  have he' : ((D.ramIdx L' : ℚ)) ≠ 0 := by
    have h := D.ramIdx_pos L'
    have h0 : D.ramIdx L' ≠ 0 := Nat.pos_iff_ne_zero.1 h
    exact_mod_cast h0
  have hcross := D.vq_baseChange L L'
  have hq : ((D.vq (D.extField L) (D.baseChange L) : ℚ)) * (D.ramIdx L' : ℚ)
      = ((D.vq (D.extField L') (D.baseChange L') : ℚ)) * (D.ramIdx L : ℚ) := by
    exact_mod_cast congrArg (fun n : ℕ => (n : ℚ)) hcross
  rw [potLocalHeight, potLocalHeight]
  field_simp
  linarith [hq]

/-! ## Proposition 3.4 —— Faltings 高さと無限遠因子 -/

/-- **[GenEll] Proposition 3.4**(Faltings Heights and the Divisor at Infinity)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

任意の `ϵ ∈ ℝ_{>0}` に対し `M_ell(ℚ̄)` 上で
`deg_∞ ≲ ht_∞ ≲ 12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht_∞`。
とくに `C ∈ ℝ` なら **`{[E] ∈ M_ell(ℚ̄)^{≤d} : ht^Falt([E]) ≤ C}` は有限**。

★**Faltings 高さが初めて出る場所**である。証明は 4 行で、
最初の `≲` は `Proposition 1.6` の証明から、残りは **[Silv2] Proposition 2.1** と
**[FC] Chapter V, Proposition 4.5**(無限遠での計量の**対数的特異性**)から。

## ★★★★★★★★ 2026-08-26 の逸脱の記録——`BDle` ではなく `BDge` で書く

以前は 3 つの `≲` を `BDle`(原典 `Definition 1.2, (ii)` の**印字どおり**)で書いていた。
★★しかし `Check/GenEll/Section3NotProvable.lean` の `bdle_chain_forces_bounded` が示すとおり、
**2 番目と 3 番目を `BDle` で読むと `ht_∞` が上に有界になってしまう**
(`(1+ε)ht_∞ ≤ ht_∞ + C` から `ε·ht_∞ ≤ C`)。
★★★モジュライ上の高さは上に有界でないので、その読みは意図されたものではない。

★★★★`Found/GenEll/BDClass.lean` の docstring がすでに定めているとおり
「**abc の主張を書くときは `BDge`(原文の印字の `≳`)を使う**」に従い、
本命題も `BDge` で書く。★これは**逸脱**であり、ここに記録する。 -/
theorem prop_3_4 (D : EllModuliData) (eps : ℝ) (heps : 0 < eps) :
    BDge D.degInf D.htInf
  ∧ BDge D.htInf (fun x => 12 * (1 + eps) * D.faltingsHeight x)
  ∧ BDge (fun x => 12 * (1 + eps) * D.faltingsHeight x) (fun x => (1 + eps) * D.htInf x)
  ∧ (∀ (C : ℝ) (d : ℕ), 0 < d → {x ∈ D.degLe d | D.faltingsHeight x ≤ C}.Finite) := by
  obtain ⟨M, hM⟩ := D.htInf_bdeq_faltings
  obtain ⟨B, hB⟩ := D.faltingsHeight_bddBelow
  refine ⟨D.degInf_le_htInf, ⟨M - 12 * eps * B, fun x => ?_⟩,
    ⟨(1 + eps) * M, fun x => ?_⟩, D.northcott⟩
  · have h1 := abs_le.1 (hM x)
    have h2 := hB x
    nlinarith [h1.1, h1.2]
  · have h1 := abs_le.1 (hM x)
    nlinarith [h1.1, h1.2]

/-! ## Lemma 3.5 —— l-捩れの大域階数 1 部分群 -/

/-- **[GenEll] Lemma 3.5**(Global Rank One Subgroups of l-Torsion)。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`l` が `E` の(悪い)乗法還元の全素点での局所高さと素なら
**`(1/(12(1+ϵ)))·l·deg_∞(E) ≤ ht^Falt(E) + 2log(l) + C`**。

★★**`C` は `ϵ` に依り得るが `E`, `F`, `H_F`, `l` に依らない**——
原文が明記している。したがって `∃ C` は `∀ E, ∀ l` の**外側**に置かねばならない。
★ここを取り違えると主張が別物になるので、量化子の順序が本 statement の要点である。 -/
theorem lemma_3_5 (D : EllModuliData) (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (E : D.Curve) (l : ℕ), Nat.Prime l →
      D.HasLCyclic E l → D.PrimeToLocalHeights E l →
      (1 / (12 * (1 + eps))) * (l : ℝ) * D.degInf (D.cls E)
        ≤ D.faltingsHeight (D.cls E) + 2 * Real.log l + C := by
  obtain ⟨_, h2, _, _⟩ := prop_3_4 D eps heps
  obtain ⟨C₁, hC₁⟩ := D.degInf_le_htInf
  obtain ⟨C₂, hC₂⟩ := h2
  obtain ⟨C₀, hC₀⟩ := D.faltingsHeight_quotLCyclic
  have hpos : (0:ℝ) < 12 * (1 + eps) := by nlinarith
  refine ⟨C₀ + (C₁ + C₂) / (12 * (1 + eps)), fun E l hl hcyc hprime => ?_⟩
  set E' := D.quotLCyclic E l with hE'
  -- deg∞(E′) = l·deg∞(E)
  have hdeg : D.degInf (D.cls E') = (l : ℝ) * D.degInf (D.cls E) :=
    D.degInf_quotLCyclic E l hl hcyc hprime
  -- deg∞(E′) ≤ 12(1+ε)·ht^Falt(E′) + (C₁ + C₂)
  have hchain : D.degInf (D.cls E')
      ≤ 12 * (1 + eps) * D.faltingsHeight (D.cls E') + (C₁ + C₂) := by
    have ha := hC₁ (D.cls E')
    have hb := hC₂ (D.cls E')
    linarith
  -- ht^Falt(E′) ≤ ht^Falt(E) + 2log(l) + C₀
  have hfal := hC₀ E l hl hcyc
  rw [hdeg] at hchain
  rw [mul_assoc, one_div_mul_eq_div, div_le_iff₀ hpos]
  have hexp : (D.faltingsHeight (D.cls E) + 2 * Real.log l
        + (C₀ + (C₁ + C₂) / (12 * (1 + eps)))) * (12 * (1 + eps))
      = 12 * (1 + eps) * (D.faltingsHeight (D.cls E) + 2 * Real.log l + C₀) + (C₁ + C₂) := by
    field_simp
    ring
  rw [hexp]
  nlinarith [hchain, hfal]

/-! ## Lemma 3.6 —— 初等的な評価 -/

/-- **[GenEll] Lemma 3.6**(An Elementary Estimate)。

原文 (GenEll p.17):
> Lemma 3.6. (An Elementary Estimate) Let

`ϵ ∈ ℝ_{>0}` に対し**定数 `C₀ ∈ ℝ_{>0}` が存在して**、
`y ≥ 1` かつ `x ≥ C₀·y^{1+ϵ}` なる全ての `x, y ∈ ℝ` について **`x ≥ y·log(x)`**。

★★**これは純粋な実解析であり、Arakelov も Galois 表現も Tate 曲線も要らない。**
原文の証明も 1 行——「`x^{1/(1+ϵ)}·log(x)/x = log(x)·x^{−ϵ/(1+ϵ)} → 0`(`x → ∞`)という
よく知られた初等的事実から直ちに従う」(p.18 目視確認 2026-08-16)。

★★**本 statement は `sorry` ではない**——`Found/GenEll/Elementary.lean` の
実装をそのまま参照している。§3 の 9 項目のうち**実装まで済んでいる唯一の項目**である。 -/
theorem lemma_3_6 (eps : ℝ) (heps : 0 < eps) :
    ∃ C₀ : ℝ, 0 < C₀ ∧ ∀ x y : ℝ, 1 ≤ y → C₀ * y ^ (1 + eps) ≤ x → y * Real.log x ≤ x :=
  ABC3.Found.GenEll.lemma_3_6 eps heps

/-! ## Lemma 3.7 —— 有限例外集合 -/

set_option maxHeartbeats 1600000 in
/-- **[GenEll] Lemma 3.7**(Finite Exceptional Sets)。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

`K_V` を compactly bounded subset、`ϵ > 0` とすると、
**定数 `C > 0` と Galois-finite な `Exc` が存在して**次を満たす:
`E_L` が半安定還元を持ち `d ≝ [L:ℚ]`、`l` 素数として、
(a) `l ≥ 100d·(ht^Falt([E_L]) + C·d^ϵ)` かつ `E_L` が乗法還元の素点を 1 つ以上持つ、
(b) `[E_L] ∈ K_V` かつ `l` が乗法還元の全素点での局所高さと素、
のいずれかを仮定すると:
(a) なら `l` は局所高さ**より大きい**。
(b) かつ `[E_L] ∉ Exc` なら `E_L` は乗法還元の素点を 1 つ以上持つ。
(a) か (b) が成り立ち、さらに `E_L` が l-cyclic 部分群スキームを持つなら **`[E_L] ∈ Exc`**。

★`Theorem 3.8` の条件 (a)(b) は**この補題の (a)(b) をそのまま引き継いでいる**
(係数だけが `100d` → `23040·100d`)。

★★**`HasMultRed`(乗法還元)と `HasPotMultRed`(潜在的乗法還元)を使い分けている**——
本補題は前者、`Theorem 3.8` は後者。`Interface` で別のフィールドにしてあるのはこのためである。 -/
theorem lemma_3_7 (D : EllModuliData) (KV : Set D.EllClass) (hKV : D.CompactlyBounded KV)
    (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, 0 < C ∧ ∃ Exc : Set D.EllClass, D.GaloisFinite Exc ∧
      ∀ (E : D.Curve) (l : ℕ), Nat.Prime l → D.SemiStable E →
        ∀ (condA condB : Prop),
          (condA ↔ (100 * (D.degOfDefinition E : ℝ)
                      * (D.faltingsHeight (D.cls E) + C * (D.degOfDefinition E : ℝ) ^ eps)
                        ≤ (l : ℝ) ∧ D.HasMultRed E)) →
          (condB ↔ (D.cls E ∈ KV ∧ D.PrimeToLocalHeights E l)) →
          (condA → D.PrimeToLocalHeights E l)
        ∧ (condB → D.cls E ∉ Exc → D.HasMultRed E)
        ∧ ((condA ∨ condB) → D.HasLCyclic E l → D.cls E ∈ Exc) := by
  -- `Proposition 3.4` を `ε₀ = 1/6`(すなわち `12(1+ε₀) = 14`)で使う
  obtain ⟨h1, h2, _, _⟩ := prop_3_4 D (1/6) (by norm_num)
  obtain ⟨C₁, hC₁⟩ := h1
  obtain ⟨C₂, hC₂⟩ := h2
  obtain ⟨B, hB⟩ := D.faltingsHeight_bddBelow
  set A : ℝ := C₁ + C₂ with hAdef
  have hdeg : ∀ x, D.degInf x ≤ 14 * D.faltingsHeight x + A := by
    intro x
    have ha := hC₁ x
    have hb := hC₂ x
    norm_num at hb
    linarith
  have hLlo : (0.69 : ℝ) ≤ Real.log 2 := by linarith [Real.log_two_gt_d9]
  have hLhi : Real.log 2 ≤ (0.70 : ℝ) := by linarith [Real.log_two_lt_d9]
  refine ⟨|A| + 100 * |B| + 1, by positivity,
    D.lcyclicExc ∪ D.noMultRedExc KV,
    D.galoisFinite_union _ _ D.galoisFinite_lcyclicExc (D.galoisFinite_noMultRedExc KV hKV), ?_⟩
  intro E l hl hss condA condB hcA hcB
  have hAimp : condA → D.PrimeToLocalHeights E l := by
    intro hc
    rw [hcA] at hc
    obtain ⟨hle, _⟩ := hc
    refine D.primeToLocalHeights_of_lt E l hl hss ?_
    set d : ℝ := (D.degOfDefinition E : ℝ) with hddef
    set F : ℝ := D.faltingsHeight (D.cls E) with hFdef
    set C : ℝ := |A| + 100 * |B| + 1 with hCdef
    set P : ℝ := d ^ eps with hPdef
    set L : ℝ := Real.log 2 with hLdef
    have hd1 : (1 : ℝ) ≤ d := by
      rw [hddef]
      exact_mod_cast D.degOfDefinition_pos E
    have hdpos : (0 : ℝ) < d := by linarith
    have hP1 : (1 : ℝ) ≤ P := by
      rw [hPdef]
      exact Real.one_le_rpow hd1 heps.le
    have hCpos : (0 : ℝ) < C := by
      rw [hCdef]
      positivity
    have hAle : A ≤ |A| := le_abs_self A
    have hAnn : (0 : ℝ) ≤ |A| := abs_nonneg A
    have hBnn : (0 : ℝ) ≤ |B| := abs_nonneg B
    have hBge : -|B| ≤ B := neg_abs_le B
    have hGlow : d * B ≤ d * F := by nlinarith [hB (D.cls E)]
    have hkey1 : d * D.degInf (D.cls E) ≤ 14 * (d * F) + d * A := by
      nlinarith [hdeg (D.cls E)]
    have hdP : d ≤ d * P := by nlinarith
    have hdPnn : (0 : ℝ) ≤ d * P := by nlinarith
    have hkey3 : 14 * (d * F) + d * A < 100 * d * (F + C * P) * L := by
      have hexp : 100 * d * (F + C * P) * L
          = 100 * L * (d * F) + 100 * L * C * (d * P) := by ring
      rw [hexp]
      have e2 : d * A ≤ d * |A| := by nlinarith
      have hX : (0 : ℝ) ≤ C * (d * P) := by nlinarith
      have e4a : 69 * (C * d) ≤ 69 * (C * (d * P)) := by nlinarith
      have e4b : 69 * (C * (d * P)) ≤ 100 * L * (C * (d * P)) := by nlinarith [hLlo, hX]
      have e4 : 69 * (C * d) ≤ 100 * L * C * (d * P) := by linarith
      rcases le_or_gt 0 (d * F) with hG | hG
      · have e1 : 14 * (d * F) ≤ 100 * L * (d * F) := by nlinarith
        have e3 : d * |A| < 69 * (C * d) := by
          rw [hCdef]
          nlinarith
        linarith
      · have f0 : -(d * F) ≤ d * |B| := by nlinarith
        have f1 : 14 * (d * F) - 100 * L * (d * F) ≤ 56 * (d * |B|) := by nlinarith
        have f3 : 56 * (d * |B|) + d * |A| < 69 * (C * d) := by
          rw [hCdef]
          nlinarith
        linarith
    have hL0 : (0 : ℝ) ≤ L := by linarith
    have hkey2 : 100 * d * (F + C * P) * L ≤ (l : ℝ) * L := by nlinarith
    linarith
  refine ⟨hAimp, ?_, ?_⟩
  · intro hc hnot
    rw [hcB] at hc
    by_contra hmult
    exact hnot (Or.inr (D.mem_noMultRedExc KV E hc.1 hmult))
  · intro hor hcyc
    have hpr : D.PrimeToLocalHeights E l := by
      rcases hor with h | h
      · exact hAimp h
      · rw [hcB] at h
        exact h.2
    exact Or.inl (D.mem_lcyclicExc E l hl hss hcyc hpr)

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def lemma_3_1.src : Source :=
  { paper := "GenEll", pdfPage := 13, item := "Lemma 3.1",
    sectionId := "genell-lemma-3-1" }

/-- ★(i)(ii)(iii) は実装済みなので、残るのは (iv) だけである。 -/
def lemma_3_1.needs : List ProofObligation :=
  [ .citation "[Serre]" "Abelian l-adic Representations and Elliptic Curves, Chapter IV, §3.4, Lemma 3(pro-l 群の Frattini 型の議論)"
      (.absent "0_Source に [Serre] は無く、mathlib にも pro-l 群の対応する補題は無い(2026-08-16、`IsProPGroup|pro-p` で 2 件ヒットするがいずれも別物)") 14,
    .implicitStep
      "★(i)(ii)(iii) は Found/GenEll/Lemma31.lean に実装済み(sorry 無し)。本 statement が sorry なのは (iv) だけによる" 14,
    .implicitStep
      "★(iv) の『閉部分群』は ℤ_l の位相から GL₂(ℤ_l) に入る位相での閉性である。**落としてはいない**が、位相群としての性質(コンパクト性・profinite 性)は使っていない" 14 ]

def lemma_3_2.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Lemma 3.2",
    sectionId := "genell-lemma-3-2" }

/-- ★★★★★★**2026-08-26 の内訳**——(i) は `Interface` の
`stableLine_dvd_or_cyclotomic` を**そのまま受けている**(原文が証明を与えず
[FC] Ch. III, Cor. 7.3 に帰しているため)。
★一方 (ii) は **`q_{E′} = q_E^l`(`vq_quotMu`)と `deg_∞ = v_K(q_E)·log #(O_K/𝔪)`
(`degInf_eq`)から導出している**——原文の『hence』の中身である。 -/
def lemma_3_2.needs : List ProofObligation :=
  [ .citation "[FC]" "Degenerations of Abelian Varieties, Chapter III, Corollary 7.3(完全列 0 → 𝔽_l(1) → M_l(E) → 𝔽_l → 0)"
      (.absent "mathlib に Tate 曲線・Tate twist・M_l(E) はいずれも無い(2026-08-16、EllipticCurve/ 配下の全宣言名を確認)") 15,
    .implicitStep
      "原文は『we have the following well-known result』とだけ書き、証明を与えていない。★(i) は「完全列の拡大類が q_E の l 乗根の抽出で得られる」ことから従うと読めるが、その段は書かれていない" 15,
    .implicitStep
      "★statement の語彙(Tate 母数 q_E・半安定還元・M_l(E)・μ_l)を Interface/GenEll/TateLocal.lean に posit した。**我々は作っていない**" 15 ]

def potLocalHeight_indep.src : Source :=
  { paper := "GenEll", pdfPage := 16, item := "Remark 3.3.1",
    sectionId := "genell-rem-3-3-1" }

def potLocalHeight_indep.needs : List ProofObligation :=
  [ .implicitStep
      "原文は『one verifies immediately that this definition is independent of the choice of L』の 1 文。★畳んでいるのは局所高さが有限拡大で分岐指数倍になること(v_L(q_E) = e(L/K)·v_K(q_E))である" 16,
    .implicitStep
      "★★★★★★**2026-08-26 に閉じた**。その実質は `Found/GenEll/LocalHeightRamified.lean` の `ordAt_liesOver`(mathlib の `HeightOneSpectrum.valuation_liesOver` から導いた)であり、`Interface/GenEll/TateLocal.lean` の `vq_baseChange` として仕様に出してある。★本定理はそこからの初等的な導出(交差乗法)である" 16 ]

def localHeight.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3",
    sectionId := "genell-def-3-3" }

def localHeight_pos.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3",
    sectionId := "genell-def-3-3" }

/-- ★**空リストは省略ではなく主張である**——正であることは `Interface` が
`vq_pos` として持っているので、原文の何にも依拠しない。

★ただし **`vq_pos` を posit したこと自体が仕事の先送りである**——
原文は「the **positive** integer `v_K(q_E)`」と書いており、
その正性は Tate 母数が `𝔪_K` に属することから来る。
`Interface/GenEll/TateLocal.lean` の `waiting` がそれを名指ししている。 -/
def localHeight_pos.needs : List ProofObligation := []

def potLocalHeight.src : Source :=
  { paper := "GenEll", pdfPage := 16, item := "Remark 3.3.1",
    sectionId := "genell-rem-3-3-1" }

def prop_3_4.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Proposition 3.4",
    sectionId := "genell-prop-3-4" }

/-- ★★★★★★**2026-08-26 に閉じた**——`Interface/GenEll/EllModuli.lean` に
原文が実際に引く 4 つ(`degInf_le_htInf`・`htInf_bdeq_faltings`・
`faltingsHeight_bddBelow`・`northcott`)を欄に出し、
**`ε` の入った 2 つの不等式はそこから初等的に導いた**。 -/
def prop_3_4.needs : List ProofObligation :=
  [ .citation "[Silv2]" "Proposition 2.1(ht_∞ と ht^Falt の比較)"
      (.absent "mathlib に Faltings 高さは無い(`Arakelov` 0 件、`arithmetic line bundle` 0 件、2026-08-16 実測)") 17,
    .citation "[FC]" "Degenerations of Abelian Varieties, Chapter V, Proposition 4.5(無限遠での計量の対数的特異性)"
      (.absent "mathlib に complex analytic space が 0 件なので、計量の特異性を述べる場所自体が無い(2026-08-16 実測)") 17,
    .citation "[Silv1]" "Proposition 8.2"
      (.unmeasured) 17,
    .otherPaper "[GenEll]" "Proposition 1.6(の証明——最初の ≲ がそこから従う)" 9,
    .otherPaper "[GenEll]" "Proposition 1.4, (iv)(Northcott の有限性)" 6 ]

def lemma_3_5.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Lemma 3.5",
    sectionId := "genell-lemma-3-5" }

/-- ★★★★★★**2026-08-26 に閉じた**——`Interface/GenEll/EllModuli.lean` に
`quotLCyclic`・`degInf_quotLCyclic`(`Lemma 3.2` の大域版)・
`faltingsHeight_quotLCyclic`([FC] Ch. I, Prop 2.7 + 複素解析の段)を欄に出し、
**`Proposition 3.4` と合わせて初等的に導いた**。 -/
def lemma_3_5.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Lemma 3.2, (i)(ii)(局所階数 1 部分群と deg_∞(E′) = l·deg_∞(E))" 15,
    .otherPaper "[GenEll]" "Proposition 3.4(Faltings 高さと無限遠因子)" 17,
    .citation "[FC]" "Degenerations of Abelian Varieties, Chapter I, Proposition 2.7(次数 l の被覆射の延長)"
      (.absent "mathlib に半安定 scheme への射の延長定理は無い(2026-08-16 実測)") 17,
    .implicitStep
      "★『(1,1)-形式を E_v 上で積分するのと (E_H)_v 上で積分するのとで l 倍だけ違う』という複素解析の段。原文は 1 文で済ませている" 17 ]

def lemma_3_6.src : Source :=
  { paper := "GenEll", pdfPage := 17, item := "Lemma 3.6",
    sectionId := "genell-lemma-3-6" }

/-- ★**空リストは省略ではなく主張である**——原文の証明は外部依存を持たない。

原文 (GenEll p.18):
> This follows immediately from the well-known elementary fact that x

★続きは `x^{1/(1+ϵ)}·log(x)/x = log(x)·x^{−ϵ/(1+ϵ)} → 0`(`x → ∞`)。
★**逐語として引けるのはここまでである**——`pdftotext` は `ϵ` を落とすので、
指数を含む部分は照合できない(2026-08-16 実測)。指数は 260 dpi 目視から読んだ。

★ここで使う「よく知られた初等的事実」は **`Real.tendsto_log_div_rpow_nhds_zero` 相当**であり、
mathlib にある。外部文献も別論文の項目も要らない。 -/
def lemma_3_6.needs : List ProofObligation := []

def lemma_3_7.src : Source :=
  { paper := "GenEll", pdfPage := 18, item := "Lemma 3.7",
    sectionId := "genell-lemma-3-7" }

/-- ★★★★★★★★**2026-08-26 に閉じた**——内訳は次のとおり。

★(a) は**完全な導出**である: `Proposition 3.4` を `ε₀ = 1/6`(すなわち `12(1+ε₀) = 14`)で
使って `deg_∞ ≤ 14·ht^Falt + A` を得、`ht^Falt ≥ B`(下に有界)と
`100·d·(ht^Falt + C·d^ε) ≤ l` から `d·deg_∞ < l·log 2` を出し、
`primeToLocalHeights_of_lt` に渡す。★★定数は `C := |A| + 100|B| + 1` と取ればよい
——`100·log 2 > 69` と `d^ε ≥ 1` が効く。

★(b)(c) は例外集合を `lcyclicExc ∪ noMultRedExc KV` と**取る**ことで導く。
★★★(c) では **(a) と (b) のどちらでも `PrimeToLocalHeights` が得られる**のが鍵である。 -/
def lemma_3_7.needs : List ProofObligation :=
  [ .otherPaper "[GenEll]" "Lemma 3.5(大域階数 1 部分群の不等式)" 17,
    .otherPaper "[GenEll]" "Proposition 3.4(『12(1+ϵ)』を 14 に取る箇所)" 17,
    .otherPaper "[GenEll]" "Proposition 1.4, (iv)(有限例外集合 Exc_d の有限性)" 6,
    .otherPaper "[GenEll]" "Example 1.3, (i)(Galois-finite)" 5,
    .implicitStep
      "★原文は『if v is any local height of E_L, then d·deg_∞([E_L]) ≥ v·log(2)』を証明なしで使う。これは局所高さの定義(v_K(q_E))と deg_∞ = log(#(O_K/(q_E))) の関係であり、素点の剰余体の位数が 2 以上であることから従うはずだが、その段は書かれていない" 18 ]

end ABC3.Skeleton.GenEll
