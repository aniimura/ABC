import ABC3.Found.PGC.HerbrandComposition

/-!
# 上付き番号付けの分岐群(Yoshida 2008 Definition 6.12 と Corollary 6.13)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Definition 6.12** と
**Corollary 6.13**(ともに物理 p.17)。構造化済み原文は
`ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/section-6.html` の
`id="def-6-12"`(`data-pdf-page="17"`, `data-item="Definition 6.12"`)と
`id="cor-6-13"`(`data-pdf-page="17"`, `data-item="Corollary 6.13"`)。

原文 (Yoshida08 p.17):
> Definition 6.12. For m ∈ R[bb]_≥0, set G^m := G_φ^−1_G(m) (the upper numbering).

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0. (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified. If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}. (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

設定は §6.2 (`#setup-6-2`):
> 6.2. The Hasse-Arf theorem. Let G = Gal(K′/K) with K′/K totally ramified as before,
> and let G ▷ H with G/H = Gal(K′′/K). For σ ∈ G, let σ[bar] = σH ∈ G/H be its image.

原典の Proof(`.txt` 1132–1143 行):
> Proof. (i): By Proposition 6.9 and Lemma 6.10(ii), we compute G^mH/H = G_{φ^{−1}_G(m)}H/H =
> (G/H)_{φ_H(φ^{−1}_G(m))} = (G/H)_{φ^{−1}_{G/H}(m)} = (G/H)^m. (ii): If G = Gal(K′K′′/K) and
> G/H = Gal(K′′/K), then G^mH/H = (G/H)^m = {id} shows G^m ⊂ H = Gal(K′K′′/K′′). Similarly
> G^m ⊂ Gal(K′K′′/K′), hence G^m = {id}. (iii): If n − 1 < φ^{−1}_G(m) ≤ n for n ∈ Z≥0, then
> G^m = G_n. …

## この定義の「本当の中身」——`φ_G` の逆写像が存在すること

原典は `φ^{−1}_G` を断りなく使うが、その存在(＝`φ_G` の全単射性)は明示されていない。
本ファイルの主題はそこである。段取りは 3 段:

1. `φ_G` は狭義単調増加(在庫 `strictMono_phiOf`。`i(1) = ⊤` の項だけから来る)。
2. `φ_G` は連続(§2 `continuous_phiOf`)。`φ_G(n) = −1 + (Σ_τ min{i(τ), n+1})/|G|` は
   有限個の `min` の和なので、`Continuous.min` と `continuous_finsetSum` の合成で出る。
3. `φ_G` は非有界(上下とも)。
   * 上: `τ = 1` の項が `min{⊤, n+1} = n+1` なので `φ_G(n) ≥ −1 + (n+1)/|G|`(`le_phiOf`)。
   * 下: ★**`i(τ) ≥ 1` なので `n ≤ 0` では全部の項が `n+1` になり、`φ_G(n) = n`**
     (`phiOf_of_nonpos`)。すなわち `φ_G` は `(−∞, 0]` の上で恒等写像である。

3 の「下」が本ファイルで一番安かった観察である。原典は `ℝ≥0` の上でしか `φ_G` を考えないが、
Y11 の `ramificationGroupReal` が `ℝ` 全体で定義されているのに合わせて `φ_G` も `ℝ` 全体で
見ると、負の側では恒等写像になる。おかげで **`φ_G : ℝ → ℝ` が全単射**になり、
中間値の定理 1 回(`intermediate_value_Icc`)で全射性が出る。
区間 `[0,∞)` への制限や `SurjOn` を扱う必要が無い。

## 逸脱の記録

1. **`m : ℝ` 全体**で `G^m` を定義した(原典は `m ∈ ℝ≥0`)。Y11 の `ramificationGroupReal`
   が `ℝ` 全体で定義されているのに合わせた。上で見たように `φ_G` は `ℝ → ℝ` の全単射
   なので `ψ_G := φ_G^{−1}` は `ℝ` 全体で定義でき、`m ≤ 0` では `G^m = G_m = ⊤` になる
   (`upperRamificationGroup_zero` はその `m = 0` の場合)。`m ≥ 0` に制限した主張は
   すべてこの一般形の特殊化である。
2. **Corollary 6.13 (i) を商群 `G ⧸ H` を作らずに述べた。** Y12 の
   `herbrand_coe_mul_coe_eq` と同じ流儀で、`G` の中の集合の等式
   `(G^m : Set G) * (H : Set G) = ((G/H)^m の引き戻し : Set G)` として述べている。
   右辺 `(G/H)^m` は `upperRamificationGroup G ϖ m`(`ϖ` は `C` の一意化元)であり、
   これが本当に商群の上の上付き分岐群であることは Y12 の
   `herbrandPhiGroup_eq_phiOf_quotient` が保証している。
   ★このため **`H` の正規性は仮定していない**(Y12 の「逸脱の記録 (3)」を引き継ぐ)。
3. **Corollary 6.13 (ii) の体論の部分を仮定に置いた。** 原典は
   「`K′K′′/K` が totally ramified」から Galois 対応で
   `Gal(K′K′′/K′) ∩ Gal(K′K′′/K′′) = {id}` を読む。本ファイルは環論の設定なので、
   これを `H ⊓ H′ = ⊥` という仮定として渡す。★ただしその仮定は
   `inf_eq_bot_of_closure_eq_top` で「`B` が 2 つの部分環 `C`, `C′` の像で生成される
   (＝合成体である)」という環論の条件から出せるようにしてある。
4. **Corollary 6.13 (iii) は本ファイルに入っていない。** Theorem 6.11(Hasse-Arf)の
   一般の可換 `G` に対する形と剰余体の位数 `q` が要り、前者は本ファイルの時点で
   木に無い(`HasseArf.lean` は巡回群・馴分岐商の場合のみ)。★新ノードとして報告した。
   その入口である「`n − 1 < ψ_G(m) ≤ n` なら `G^m = G_n`」だけは
   `upperRamificationGroup_eq_lowerRamificationGroup` として本ファイルに入れてある。

## 退化の自己検査

* `φ_G` の全単射性を **仮定に置いていない**(§4 で証明している)。逆写像は
  `Function.invFun` で作り、`φ_G ∘ ψ_G = id` は全射性から、`ψ_G ∘ φ_G = id` は
  単射性から出る。★全射性には `i(σ) ≥ 1`(`pos_ramIndex`、`huni` が要る)が要る。
  これを落とすと `φ_G` が `(−∞,0]` で恒等でなくなり、下への非有界性が崩れる。
* `m` の `≥ 0` を落としても `ψ_G(m)` は定義される(逸脱 1)。`m ≥ 0` を使うのは
  `herbrandPsiGroup_nonneg` と `exists_upperRamificationGroup_eq_lowerRamificationGroup`
  だけである。
* (i) の `G ▷ H` は **使っていない**(逸脱 2)。商群を作らないので正規性が要らない。
* (iii) の `G` 可換は本ファイルでは出番が無い(iii 自体が入っていないため)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing
open scoped Pointwise

/-! ## §1 抽象核 —— 純解析

★群も分岐も付値も出てこない。`f : ℝ → ℝ` だけの話。 -/

/-- ★★**抽象核** —— 狭義単調・連続で上にも下にも非有界な `f : ℝ → ℝ` は全射。

★段取り: `f a ≤ y ≤ f b` を取ると狭義単調性から `a ≤ b`、あとは中間値の定理
`intermediate_value_Icc` 1 回。★`Ici 0` への制限も `SurjOn` も要らない
(`φ_G` が `ℝ` 全体で全単射になるように定義してあるため)。 -/
theorem surjective_of_strictMono_continuous {f : ℝ → ℝ} (hmono : StrictMono f)
    (hcont : Continuous f) (hup : ∀ y : ℝ, ∃ b : ℝ, y ≤ f b)
    (hdown : ∀ y : ℝ, ∃ a : ℝ, f a ≤ y) : Function.Surjective f := by
  intro y
  obtain ⟨a, ha⟩ := hdown y
  obtain ⟨b, hb⟩ := hup y
  have hab : a ≤ b := hmono.le_iff_le.1 (ha.trans hb)
  obtain ⟨x, _, hx⟩ := intermediate_value_Icc hab hcont.continuousOn ⟨ha, hb⟩
  exact ⟨x, hx⟩

/-- ★**抽象核** —— 狭義単調な全射の逆写像もまた狭義単調。 -/
theorem strictMono_invFun {f : ℝ → ℝ} (hmono : StrictMono f) (hsurj : Function.Surjective f) :
    StrictMono (Function.invFun f) := by
  intro m m' hmm
  have h1 : f (Function.invFun f m) = m := Function.invFun_eq (hsurj m)
  have h2 : f (Function.invFun f m') = m' := Function.invFun_eq (hsurj m')
  exact hmono.lt_iff_lt.1 (by rw [h1, h2]; exact hmm)

/-! ## §2 抽象核 —— `phiOf` の解析的性質

★Y12 の `phiOf : (S → ℕ∞) → ℝ → ℝ`(有限型 `S` と `f : S → ℕ∞` だけの定義)に対して、
連続性・狭義単調性・非有界性・全射性をここで一気に片づける。
分岐・付値・Galois の語彙は 1 つも出てこない。 -/

/-- `x ≥ 1` かつ `r ≤ 1` なら `min{x, r} = r`。★`ℕ∞` の側の `x` は `⊤` かもしれないので
`truncENat` の場合分けを 1 回だけここで済ませる。 -/
theorem truncENat_of_le_one {x : ℕ∞} (hx : 0 < x) {r : ℝ} (hr : r ≤ 1) :
    truncENat x r = r := by
  rcases eq_or_ne x ⊤ with h | h
  · rw [h, truncENat_top]
  · rw [← ENat.coe_toNat h, truncENat_coe]
    have h1 : 1 ≤ x.toNat := by
      have h2 : (1 : ℕ∞) ≤ x := Order.one_le_iff_pos.2 hx
      simpa using ENat.toNat_le_toNat h2 h
    exact min_eq_right (le_trans hr (by exact_mod_cast h1))

/-- `r ≥ 0` なら `min{x, r} ≥ 0`。★`x` の正値性は要らない(`(x.toNat : ℝ) ≥ 0` だから)。 -/
theorem truncENat_nonneg (x : ℕ∞) {r : ℝ} (hr : 0 ≤ r) : 0 ≤ truncENat x r := by
  unfold truncENat
  split
  · exact hr
  · exact le_min (by positivity) hr

/-- ★★★**`φ` は `(−∞, 0]` の上で恒等写像**(`f > 0` が要る)。

★これが本ファイルで一番安い観察である。`n ≤ 0` では `n + 1 ≤ 1 ≤ f(τ)` なので
すべての項が `min{f(τ), n+1} = n+1` になり、`φ(n) = −1 + |S|(n+1)/|S| = n`。
★系として `φ` は下に非有界(`φ(min y 0) = min y 0 ≤ y`)。
★`n = 0` の場合が在庫の `phiOf_zero`(Lemma 6.10 (i) の前半)である。 -/
theorem phiOf_of_nonpos {S : Type*} [Fintype S] [Nonempty S] {f : S → ℕ∞}
    (hf : ∀ τ, 0 < f τ) {n : ℝ} (hn : n ≤ 0) : phiOf f n = n := by
  have hcard : (Nat.card S : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  have hsum : ∑ τ : S, truncENat (f τ) (n + 1) = (Nat.card S : ℝ) * (n + 1) := by
    rw [Finset.sum_congr rfl fun τ _ => truncENat_of_le_one (hf τ) (by linarith : n + 1 ≤ 1)]
    rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, Nat.card_eq_fintype_card]
  rw [phiOf, hsum]
  field_simp
  ring

/-- ★★**`φ` の下からの評価** `φ(n) ≥ −1 + (n+1)/|S|`(`n ≥ 0`)。

★根拠は `f(τ₀) = ⊤` の項 1 つだけ(`min{⊤, n+1} = n+1`)。他の項は非負であればよい。
★系として `φ` は上に非有界。具体層では `τ₀ = 1`(`ramIndex_one`)。 -/
theorem le_phiOf {S : Type*} [Fintype S] {f : S → ℕ∞} {τ₀ : S} (hτ₀ : f τ₀ = ⊤)
    {n : ℝ} (hn : 0 ≤ n) : -1 + (n + 1) / (Nat.card S : ℝ) ≤ phiOf f n := by
  have hpos : (0 : ℝ) < (Nat.card S : ℝ) := by
    have : Nonempty S := ⟨τ₀⟩
    exact_mod_cast Nat.card_pos
  have hle : n + 1 ≤ ∑ τ : S, truncENat (f τ) (n + 1) := by
    have := Finset.single_le_sum (f := fun τ : S => truncENat (f τ) (n + 1))
      (fun τ _ => truncENat_nonneg (f τ) (by linarith)) (Finset.mem_univ τ₀)
    rwa [hτ₀, truncENat_top] at this
  simp only [phiOf, add_le_add_iff_left]
  gcongr

/-- ★**`r ↦ min{x, r}` は連続**。`x = ⊤` なら恒等、そうでなければ `min` の片側が定数。 -/
theorem continuous_truncENat (x : ℕ∞) : Continuous (truncENat x) := by
  unfold truncENat
  by_cases hx : x = ⊤
  · simp only [hx, if_true]; exact continuous_id
  · simp only [hx, if_false]; exact continuous_const.min continuous_id

/-- ★★★**`φ` は連続** —— 有限個の `min` の和を定数で割ったものだから。

★原典は Lemma 6.10 (ii) の証明で「`φ` is continuous and piecewise linear」と
一言で済ませている。区分線形性は使わず、連続性だけをここで取る。 -/
theorem continuous_phiOf {S : Type*} [Fintype S] (f : S → ℕ∞) : Continuous (phiOf f) := by
  unfold phiOf
  exact continuous_const.add
    (Continuous.div_const
      (continuous_finsetSum _ fun i _ => (continuous_truncENat (f i)).comp
        (continuous_id.add continuous_const)) _)

/-- ★★**`φ` は狭義単調増加**(`f(τ₀) = ⊤` なる `τ₀` が 1 つあればよい)。

★在庫の `strictMono_herbrandPhi`(部分群版)の `phiOf` 版。狭義性は `τ₀` の項
`min{⊤, n+1} = n+1` だけから来る。 -/
theorem strictMono_phiOf {S : Type*} [Fintype S] {f : S → ℕ∞} (hf : ∃ τ₀, f τ₀ = ⊤) :
    StrictMono (phiOf f) := by
  obtain ⟨τ₀, hτ₀⟩ := hf
  have hpos : (0 : ℝ) < (Nat.card S : ℝ) := by
    have : Nonempty S := ⟨τ₀⟩
    exact_mod_cast Nat.card_pos
  intro a b hab
  simp only [phiOf, add_lt_add_iff_left]
  refine (div_lt_div_iff_of_pos_right hpos).mpr ?_
  refine Finset.sum_lt_sum (fun τ _ => monotone_truncENat _ (by linarith))
    ⟨τ₀, Finset.mem_univ _, ?_⟩
  rw [hτ₀, truncENat_top, truncENat_top]
  linarith

/-- ★★★★**抽象核の到達点** —— `f > 0` かつ `f(τ₀) = ⊤` なら `phiOf f : ℝ → ℝ` は全射。

★上への非有界性は `le_phiOf`(`b := max 0 (|S|(y+1))` を取る)、
下への非有界性は `phiOf_of_nonpos`(`a := min y 0` を取る)から。 -/
theorem surjective_phiOf {S : Type*} [Fintype S] {f : S → ℕ∞} (hpos : ∀ τ, 0 < f τ)
    {τ₀ : S} (hτ₀ : f τ₀ = ⊤) : Function.Surjective (phiOf f) := by
  have hne : Nonempty S := ⟨τ₀⟩
  have hN : (0 : ℝ) < (Nat.card S : ℝ) := by exact_mod_cast Nat.card_pos (α := S)
  refine surjective_of_strictMono_continuous (strictMono_phiOf ⟨τ₀, hτ₀⟩)
    (continuous_phiOf f) (fun y => ?_) (fun y => ?_)
  · refine ⟨max 0 ((Nat.card S : ℝ) * (y + 1)), le_trans ?_ (le_phiOf hτ₀ (le_max_left _ _))⟩
    have hb : (Nat.card S : ℝ) * (y + 1) ≤ max 0 ((Nat.card S : ℝ) * (y + 1)) := le_max_right _ _
    have h2 : (y + 1) ≤ (max 0 ((Nat.card S : ℝ) * (y + 1)) + 1) / (Nat.card S : ℝ) := by
      rw [le_div_iff₀ hN]
      nlinarith
    linarith
  · exact ⟨min y 0, le_of_eq_of_le (phiOf_of_nonpos hpos (min_le_right _ _)) (min_le_left _ _)⟩

/-- `phiOf f : ℝ → ℝ` は全単射。 -/
theorem bijective_phiOf {S : Type*} [Fintype S] {f : S → ℕ∞} (hpos : ∀ τ, 0 < f τ)
    {τ₀ : S} (hτ₀ : f τ₀ = ⊤) : Function.Bijective (phiOf f) :=
  ⟨(strictMono_phiOf ⟨τ₀, hτ₀⟩).injective, surjective_phiOf hpos hτ₀⟩

/-! ## §3 抽象核 —— 純群論(Corollary 6.13 (ii) の骨)

★★分岐・付値・Galois の語彙が 1 つも出てこない。 -/

/-- ★**`S·H = H` なら `S ≤ H`**。`σ = σ · 1` を見るだけ。

★原文「G^mH/H = (G/H)^m = {id} shows G^m ⊂ H」がこれである。 -/
theorem le_of_coe_mul_coe_eq {G : Type*} [Group G] {S H : Subgroup G}
    (h : (S : Set G) * (H : Set G) = (H : Set G)) : S ≤ H := by
  intro σ hσ
  have hm : σ ∈ (S : Set G) * (H : Set G) := ⟨σ, hσ, 1, H.one_mem, mul_one σ⟩
  rw [h] at hm
  exact hm

/-- ★**2 つの部分群に入り、その交わりが自明なら自明**。

★原文「G^m ⊂ H かつ G^m ⊂ Gal(K′K′′/K′)、hence G^m = {id}」がこれである。 -/
theorem eq_bot_of_le_of_inf_eq_bot {G : Type*} [Group G] {S H H' : Subgroup G}
    (h1 : S ≤ H) (h2 : S ≤ H') (hHH : H ⊓ H' = ⊥) : S = ⊥ :=
  le_bot_iff.1 (hHH ▸ le_inf h1 h2)

/-- ★★★**原文「K′K′′/K」(合成体)の環論版** ——
`B` が 2 つの部分環 `C`, `C′` の像で生成されるなら、
`C` に自明に作用する部分群 `H` と `C′` に自明に作用する部分群 `H′` の交わりは自明。

★原典は Galois 対応で `Gal(K′K′′/K′) ∩ Gal(K′K′′/K′′) = Gal(K′K′′/K′K′′) = {id}` と読む。
本補題はそれを「`σ` の固定部分環(`RingHom.eqLocus`)が 2 つの像を含むから `⊤`」
という環論の言い方に置き換えたものである。★`FaithfulSMul G B` が要る。 -/
theorem inf_eq_bot_of_closure_eq_top {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [FaithfulSMul G B] {H H' : Subgroup G}
    {C C' : Type*} [CommRing C] [CommRing C'] [Algebra C B] [Algebra C' B]
    [MulSemiringAction G C] [MulSemiringAction G C']
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hcomp' : ∀ (ρ : G) (c : C'), algebraMap C' B (ρ • c) = ρ • algebraMap C' B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hH'triv : ∀ ρ ∈ H', ∀ c : C', ρ • c = c)
    (hgen : Subring.closure (Set.range (algebraMap C B) ∪ Set.range (algebraMap C' B)) = ⊤) :
    H ⊓ H' = ⊥ := by
  rw [eq_bot_iff]
  intro σ hσ
  obtain ⟨hσH, hσH'⟩ := Subgroup.mem_inf.1 hσ
  have hsub : Subring.closure (Set.range (algebraMap C B) ∪ Set.range (algebraMap C' B))
      ≤ RingHom.eqLocus (MulSemiringAction.toRingHom G B σ) (RingHom.id B) := by
    refine Subring.closure_le.2 ?_
    rintro x (⟨c, rfl⟩ | ⟨c, rfl⟩)
    · show σ • algebraMap C B c = algebraMap C B c
      rw [← hcomp σ c, hHtriv σ hσH c]
    · show σ • algebraMap C' B c = algebraMap C' B c
      rw [← hcomp' σ c, hH'triv σ hσH' c]
  rw [hgen, top_le_iff] at hsub
  have hall : ∀ b : B, σ • b = b := fun b => by
    have hb : b ∈ RingHom.eqLocus (MulSemiringAction.toRingHom G B σ) (RingHom.id B) := by
      rw [hsub]; trivial
    exact hb
  exact Subgroup.mem_bot.2 (eq_of_smul_eq_smul (M := G) (α := B) (by simpa using hall))

/-! ## §4 具体層 —— `φ_G : ℝ → ℝ` は全単射

★抽象核 §2 に `f := i` を代入するだけ。`i(1) = ⊤` は `ramIndex_one`、
`i(σ) ≥ 1` は `pos_ramIndex`(`huni` が要る)。 -/

section PhiGroup

variable {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
variable {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G]

/-- ★★**`φ_G` は狭義単調増加**。★仮定は何も要らない(`i(1) = ⊤` だけで出る)。 -/
theorem strictMono_herbrandPhiGroup (α : B) : StrictMono (herbrandPhiGroup G α) :=
  strictMono_phiOf ⟨1, ramIndex_one α⟩

/-- ★**`φ_G` は `(−∞, 0]` の上で恒等写像**。★`φ_G(0) = 0`(在庫 `herbrandPhiGroup_zero`)の
一般化であり、下への非有界性の根拠でもある。 -/
theorem herbrandPhiGroup_of_nonpos {α : B} (huni : maximalIdeal B = Ideal.span {α}) {n : ℝ}
    (hn : n ≤ 0) : herbrandPhiGroup G α n = n :=
  phiOf_of_nonpos (fun σ : G => pos_ramIndex huni σ) hn

/-- ★★★★**本ファイルの土台** —— `φ_G : ℝ → ℝ` は全射。

★★原典は `φ^{−1}_G` を断りなく使うが、その存在の根拠はここにしか無い。 -/
theorem surjective_herbrandPhiGroup {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    Function.Surjective (herbrandPhiGroup G α) :=
  surjective_phiOf (fun σ : G => pos_ramIndex huni σ) (ramIndex_one (G := G) α)

/-- `φ_G : ℝ → ℝ` は全単射。 -/
theorem bijective_herbrandPhiGroup {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    Function.Bijective (herbrandPhiGroup G α) :=
  ⟨(strictMono_herbrandPhiGroup α).injective, surjective_herbrandPhiGroup huni⟩

end PhiGroup

/-! ## §5 Definition 6.12 —— 上付き番号付け -/

def herbrandPsiGroup.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Definition 6.12", sectionId := "def-6-12" }

/-- ★★★★**原典の `φ^{−1}_G`** —— Herbrand 関数 `φ_G : ℝ → ℝ` の逆写像 `ψ_G`。

原文 (Yoshida08 p.17):
> Definition 6.12. For m ∈ R[bb]_≥0, set G^m := G_φ^−1_G(m) (the upper numbering).

`Function.invFun` で作ってある。逆写像であること(両側)は
`herbrandPhiGroup_herbrandPsiGroup`(全射性が要る)と
`herbrandPsiGroup_herbrandPhiGroup`(単射性だけでよい)にある。

★★`φ_G` が全単射であること自体は §4 で証明してある。仮定には置いていない。 -/
noncomputable def herbrandPsiGroup {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] (G : Type*) [Group G] [MulSemiringAction G B] [Fintype G]
    (α : B) (m : ℝ) : ℝ :=
  Function.invFun (herbrandPhiGroup G α) m

section Psi

variable {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
variable {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G]

/-- `ψ_G(φ_G(n)) = n`。★単射性(＝狭義単調性)だけでよいので仮定が要らない。 -/
theorem herbrandPsiGroup_herbrandPhiGroup (α : B) (n : ℝ) :
    herbrandPsiGroup G α (herbrandPhiGroup G α n) = n :=
  Function.leftInverse_invFun (strictMono_herbrandPhiGroup α).injective n

/-- `φ_G(ψ_G(m)) = m`。★こちらは全射性(§4)が要る。 -/
theorem herbrandPhiGroup_herbrandPsiGroup {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (m : ℝ) : herbrandPhiGroup G α (herbrandPsiGroup G α m) = m :=
  Function.invFun_eq (surjective_herbrandPhiGroup huni m)

/-- `ψ_G` も狭義単調増加。 -/
theorem strictMono_herbrandPsiGroup {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    StrictMono (herbrandPsiGroup G α) :=
  strictMono_invFun (strictMono_herbrandPhiGroup α) (surjective_herbrandPhiGroup huni)

/-- `ψ_G(0) = 0`(`φ_G(0) = 0` の裏返し)。 -/
theorem herbrandPsiGroup_zero {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    herbrandPsiGroup G α 0 = 0 := by
  have h := herbrandPsiGroup_herbrandPhiGroup (G := G) α 0
  rwa [herbrandPhiGroup_zero huni] at h

/-- ★**退化検査** `m ≥ 0 ⟹ ψ_G(m) ≥ 0`。
原典の `m ∈ ℝ≥0` が `φ^{−1}_G(m) ∈ ℝ≥0` を保証している、ということ。 -/
theorem herbrandPsiGroup_nonneg {α : B} (huni : maximalIdeal B = Ideal.span {α}) {m : ℝ}
    (hm : 0 ≤ m) : 0 ≤ herbrandPsiGroup G α m := by
  rw [← herbrandPsiGroup_zero (G := G) huni]
  exact (strictMono_herbrandPsiGroup huni).monotone hm

end Psi

def upperRamificationGroup.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Definition 6.12", sectionId := "def-6-12" }

/-- ★★★★**Yoshida 2008 Definition 6.12(上付き番号付け)** —— `G^m := G_{φ^{−1}_G(m)}`。

原文 (Yoshida08 p.17):
> Definition 6.12. For m ∈ R[bb]_≥0, set G^m := G_φ^−1_G(m) (the upper numbering).

`G_·` は Y11 の実数添字の下付き分岐群 `ramificationGroupReal`、`φ^{−1}_G` は
`herbrandPsiGroup`。★これが pGC §2 の Definition 2.3(`RamificationFiltration p`)と
同じ対象である。

★逸脱: 原典は `m ∈ ℝ≥0` だが、ここは `m : ℝ` 全体(冒頭「逸脱の記録 1」)。
`m ≤ 0` では `G^m = ⊤` になる。 -/
noncomputable def upperRamificationGroup {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] (G : Type*) [Group G] [MulSemiringAction G B] [Fintype G]
    (α : B) (m : ℝ) : Subgroup G :=
  ramificationGroupReal α (herbrandPsiGroup G α m)

section Upper

variable {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
variable {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G]

theorem upperRamificationGroup_def (G : Type*) [Group G] [MulSemiringAction G B] [Fintype G]
    (α : B) (m : ℝ) :
    upperRamificationGroup G α m = ramificationGroupReal α (herbrandPsiGroup G α m) := rfl

@[simp] theorem mem_upperRamificationGroup {α : B} {m : ℝ} {σ : G} :
    σ ∈ upperRamificationGroup G α m ↔ RealLeENat (herbrandPsiGroup G α m + 1) (ramIndex α σ) :=
  Iff.rfl

/-- ★★**上付きと下付きの辞書** `G^{φ_G(n)} = G_n`。

★★原典の 4 つ組 `G^m = G_{φ^{−1}(m)}` を `m = φ_G(n)` で読み替えた形で、
Corollary 6.13 (i) の計算はすべてこの形で進む。 -/
theorem upperRamificationGroup_herbrandPhiGroup (α : B) (n : ℝ) :
    upperRamificationGroup G α (herbrandPhiGroup G α n) = ramificationGroupReal α n := by
  rw [upperRamificationGroup_def, herbrandPsiGroup_herbrandPhiGroup]

/-- `G^m` は `m` について反単調。 -/
theorem upperRamificationGroup_antitone {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    Antitone (upperRamificationGroup G α) := fun _ _ hab =>
  ramificationGroupReal_antitone α ((strictMono_herbrandPsiGroup huni).monotone hab)

omit [Fintype G] in
/-- `G_0 = G`(実数添字版)。★`i(σ) ≥ 1` だけで出る。 -/
theorem ramificationGroupReal_zero_eq_top {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    ramificationGroupReal (G := G) α 0 = ⊤ := by
  ext σ
  simp only [mem_ramificationGroupReal, Subgroup.mem_top, iff_true]
  have h : ((1 : ℕ) : ℝ) = (0 : ℝ) + 1 := by norm_num
  rw [← h, realLeENat_natCast_iff]
  exact Order.one_le_iff_pos.2 (pos_ramIndex huni σ)

/-- ★**`G^0 = G`** —— `ψ_G(0) = 0` と `G_0 = G` から。 -/
theorem upperRamificationGroup_zero {α : B} (huni : maximalIdeal B = Ideal.span {α}) :
    upperRamificationGroup G α 0 = ⊤ := by
  rw [upperRamificationGroup_def, herbrandPsiGroup_zero huni, ramificationGroupReal_zero_eq_top huni]

end Upper

section UpperLower

variable {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
  [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G]
  [SMulCommClass G A B]

/-- ★★**Corollary 6.13 (iii) の入口** ——
原文「If n − 1 < φ^{−1}_G(m) ≤ n for n ∈ Z≥0, then G^m = G_n」。

★Y11 の `ramificationGroupReal_eq_of_mem_Ioc`(実数添字の `G_x` が
`x ∈ (n−1, n]` で `ℕ` 添字の `G_n` に一致する)に `x := ψ_G(m)` を代入しただけ。 -/
theorem upperRamificationGroup_eq_lowerRamificationGroup {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    {m : ℝ} (n : ℕ) (h1 : (n : ℝ) - 1 < herbrandPsiGroup G α m)
    (h2 : herbrandPsiGroup G α m ≤ (n : ℝ)) :
    upperRamificationGroup G α m = lowerRamificationGroup B G n := by
  rw [upperRamificationGroup_def]
  exact ramificationGroupReal_eq_of_mem_Ioc (A := A) huni hadj n h1 h2

/-- ★★**`m ≥ 0` なら `G^m` は必ず `ℕ` 添字の下付き分岐群に一致する**(`n := ⌈ψ_G(m)⌉`)。

★Corollary 6.13 (iii) の第 1 文であり、pGC §2 が `G^m` を `G_n` に翻訳して使うときの入口。 -/
theorem exists_upperRamificationGroup_eq_lowerRamificationGroup {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    {m : ℝ} (hm : 0 ≤ m) :
    ∃ n : ℕ, upperRamificationGroup G α m = lowerRamificationGroup B G n := by
  have hnn : 0 ≤ herbrandPsiGroup G α m := herbrandPsiGroup_nonneg huni hm
  refine ⟨⌈herbrandPsiGroup G α m⌉₊,
    upperRamificationGroup_eq_lowerRamificationGroup (A := A) huni hadj _ ?_ (Nat.le_ceil _)⟩
  have h := Nat.ceil_lt_add_one hnn
  linarith

end UpperLower

/-! ## §6 Corollary 6.13 (i) -/

def upperRamification_coe_mul_coe_eq.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★**Yoshida 2008 Corollary 6.13 (i)** —— `G^m H/H = (G/H)^m`。

原文 (Yoshida08 p.17):
> Corollary 6.13. (i) If G ⊲ H, then G^mH/H = (G/H)^m for all m ∈ R[bb]_≥0.

原文の証明(1 行):
> By Proposition 6.9 and Lemma 6.10(ii), we compute G^mH/H = G_{φ^{−1}_G(m)}H/H
> = (G/H)_{φ_H(φ^{−1}_G(m))} = (G/H)_{φ^{−1}_{G/H}(m)} = (G/H)^m.

★★段取りは原文そのまま。`n := ψ_G(m)` と置くと
1. Lemma 6.10 (ii)(在庫 `herbrandPhiGroup_comp`)で `φ_{G/H}(φ_H(n)) = φ_G(n) = m`。
2. `φ_{G/H}` は単射(§4 の狭義単調性)なので `ψ_{G/H}(m) = φ_H(n)`。
   ★ここは全射性を使わないので `ϖ` が `C` の一意化元であることを仮定に足す必要が無い。
3. あとは Proposition 6.9(在庫 `herbrand_coe_mul_coe_eq`)を `n` に適用するだけ。

★逸脱: 商群 `G ⧸ H` を作らず、`G` の中の集合の等式として述べている
(冒頭「逸脱の記録 2」)。右辺 `upperRamificationGroup G ϖ m` が本当に
商群の上の `(G/H)^m` の引き戻しであることは、Y12 の
`herbrandPhiGroup_eq_phiOf_quotient` と Y11 の `herbrand_coe_mul_coe_eq` の
docstring が保証している。★そのため `H` の正規性は仮定に無い。

仮定はすべて Y11 の Proposition 6.9 / Y12 の Lemma 6.10 (ii) のものを
そのまま引き継いだ。★本定理が新たに足した仮定は 1 つも無い。 -/
theorem upperRamification_coe_mul_coe_eq {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [Fintype G] [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (m : ℝ) :
    ((upperRamificationGroup G π' m : Subgroup G) : Set G) * (H : Set G)
      = ((upperRamificationGroup G ϖ m : Subgroup G) : Set G) := by
  have huni : maximalIdeal B = Ideal.span {π'} := (irreducible_iff_uniformizer π').1 hπ'
  -- 1. Lemma 6.10 (ii): `φ_{G/H}(φ_H(ψ_G(m))) = φ_G(ψ_G(m)) = m`
  have hcomp2 : herbrandPhiGroup G ϖ (herbrandPhi π' H (herbrandPsiGroup G π' m)) = m := by
    rw [← herbrandPhiGroup_comp (A := A) hcomp hHtriv hπ' hinj hfixC hres hAC hϖ hadj hfix]
    exact herbrandPhiGroup_herbrandPsiGroup huni m
  -- 2. `φ_{G/H}` の単射性で `ψ_{G/H}(m) = φ_H(ψ_G(m))`
  have hψ : herbrandPsiGroup G ϖ m = herbrandPhi π' H (herbrandPsiGroup G π' m) := by
    conv_lhs => rw [← hcomp2]
    exact herbrandPsiGroup_herbrandPhiGroup _ _
  -- 3. Proposition 6.9
  rw [upperRamificationGroup_def, upperRamificationGroup_def, hψ]
  exact herbrand_coe_mul_coe_eq (A := A) hcomp hHtriv hπ' hinj hfixC hres hAC hϖ hadj hfix _

/-! ## §7 Corollary 6.13 (ii) -/

/-- ★★**(ii) の前半** —— `(G/H)^m = {id}` なら `G^m ⊆ H`。

原文:
> If G = Gal(K′K′′/K) and G/H = Gal(K′′/K), then G^mH/H = (G/H)^m = {id} shows
> G^m ⊂ H = Gal(K′K′′/K′′).

★`(G/H)^m = {id}` は、引き戻しの言い方では `upperRamificationGroup G ϖ m = H`。
(i) を使って `G^m · H = H` にし、§3 の `le_of_coe_mul_coe_eq` を当てるだけ。 -/
theorem upperRamificationGroup_le_of_quotient_eq_bot {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B] [FaithfulSMul G B]
    {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    {m : ℝ} (htriv : upperRamificationGroup G ϖ m = H) :
    upperRamificationGroup G π' m ≤ H := by
  refine le_of_coe_mul_coe_eq ?_
  rw [upperRamification_coe_mul_coe_eq (A := A) hcomp hHtriv hπ' hinj hfixC hres hAC hϖ hadj hfix,
    htriv]

def upperRamificationGroup_eq_bot_of_two_quotients.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Corollary 6.13", sectionId := "cor-6-13" }

/-- ★★★★**Yoshida 2008 Corollary 6.13 (ii)** —— 合成体でも上付き分岐が消える。

原文 (Yoshida08 p.17):
> (ii) Let K′/K and K′′/K be two Galois extensions with K′K′′/K totally ramified.
> If Gal(K′/K)^m = Gal(K′′/K)^m = {id} for m ∈ R[bb]_≥0, then Gal(K′K′′/K)^m = {id}.

`G = Gal(K′K′′/K)`、`H = Gal(K′K′′/K′′)`、`H′ = Gal(K′K′′/K′)`。
`C`(一意化元 `ϖ`)が `O_{K′′}`、`C′`(一意化元 `ϖ′`)が `O_{K′}` に対応する。
仮定 `htriv` / `htriv′` が原文の `Gal(K′′/K)^m = Gal(K′/K)^m = {id}`(引き戻しの形)である。

★★逸脱(冒頭「逸脱の記録 3」): 原文の `K′K′′/K` が totally ramified という条件から
Galois 対応で読む `Gal(K′K′′/K′) ∩ Gal(K′K′′/K′′) = {id}` を、
仮定 `hHH : H ⊓ H′ = ⊥` として渡している。★この仮定は §3 の
`inf_eq_bot_of_closure_eq_top`(`B` が `C` と `C′` の像で生成される＝合成体である)
から出せる。

★段取りは原文どおり: (i) を 2 回使って `G^m ≤ H` と `G^m ≤ H′` を出し、
交わりが自明であることから `G^m = ⊥`。 -/
theorem upperRamificationGroup_eq_bot_of_two_quotients {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B] [FaithfulSMul G B]
    {H H' : Subgroup G} [Fintype H] [Fintype H'] {π' π'' π''' : B}
    {C C' : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    [CommRing C'] [IsDomain C'] [IsDiscreteValuationRing C'] [Algebra C' B]
    [MulSemiringAction G C']
    (hπ' : Irreducible π')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (hcomp' : ∀ (ρ : G) (c : C'), algebraMap C' B (ρ • c) = ρ • algebraMap C' B c)
    (hH'triv : ∀ ρ ∈ H', ∀ c : C', ρ • c = c)
    (hinj' : Function.Injective (algebraMap C' B))
    (hfixC' : ∀ b : B, (∀ ρ ∈ H', ρ • b = b) → ∃ c : C', algebraMap C' B c = b)
    (hres' : ∀ b : B, ∃ c : C', b - algebraMap C' B c ∈ maximalIdeal B)
    (hAC' : ∀ a : A, ∃ c : C', algebraMap C' B c = algebraMap A B a)
    {ϖ' : C'} (hϖ' : algebraMap C' B ϖ' = π''')
    (hfix' : ∀ y : B, (∀ ρ ∈ H', ρ • y = y) → y ∈ Algebra.adjoin A ({π'''} : Set B))
    (hHH : H ⊓ H' = ⊥)
    {m : ℝ} (htriv : upperRamificationGroup G ϖ m = H)
    (htriv' : upperRamificationGroup G ϖ' m = H') :
    upperRamificationGroup G π' m = ⊥ :=
  eq_bot_of_le_of_inf_eq_bot
    (upperRamificationGroup_le_of_quotient_eq_bot (A := A) hcomp hHtriv hπ' hinj hfixC hres hAC
      hϖ hadj hfix htriv)
    (upperRamificationGroup_le_of_quotient_eq_bot (A := A) hcomp' hH'triv hπ' hinj' hfixC' hres'
      hAC' hϖ' hadj hfix' htriv')
    hHH

/-! ## §8 Corollary 6.13 (iii) について

★★**本ファイルには入っていない。** 原文は

> (iii) Let G be abelian. Then |G/G^m| divides (q − 1)q^m−1 for m ∈ Z[bb]_≥0.

で、証明は
「`n − 1 < φ^{−1}_G(m) ≤ n` なら `G^m = G_n`。`1 ≤ i ≤ n` の `G_i` を見ると、
Theorem 6.11 により `G_{i−1} ≠ G_i` が起きるのは `φ_G(i−1) ∈ ℤ` のときだけで、
`0 ≤ φ_G(i−1) ≤ φ_G(n−1) < m` だから `i > 1` では高々 `m − 1` 回。
Proposition 6.2 により `|G_{i−1}/G_i|` は `i = 1` のとき `q − 1` を、`i > 1` のとき `q` を割る」
というものである。要るのは

1. **Theorem 6.11(Hasse-Arf)の一般の可換 `G` に対する形** ——
   木にある `HasseArf.lean` は巡回群と馴分岐商の場合(`exists_natCast_herbrandPhiGroup_of_isCyclic`
   / `..._of_tame_quotient`)までで、Proposition 6.2 の直和分解に沿った帰納が要る。
2. **剰余体の位数 `q`** と `|G_{i−1}/G_i| ∣ q`(Proposition 6.2)。
3. 商の位数の望遠鏡積 `|G/G_n| = ∏_{i=1}^{n} |G_{i−1}/G_i|` と
   「跳びが高々 `m − 1` 回」の数え上げ。

★入口(上の第 1 文)は `upperRamificationGroup_eq_lowerRamificationGroup` /
`exists_upperRamificationGroup_eq_lowerRamificationGroup` として本ファイルに入れてある。
残りは新ノードとして報告した。 -/

end ABC3.Found.PGC
