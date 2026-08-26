import ABC3.Found.GaloisRep.WeierstrassQuotient

/-!
# Galois (G6) 第 246 ブロック —— **★★★★★★★共線性の行列式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★経路の切り替え —— 商ではなく行列式で

§9-557・§9-558 では加法定理を**商の形**

    F(z) = ℘(z+w) + ℘(z) + ℘(w) − (1/4)((℘'(z)−℘'(w))/(℘(z)−℘(w)))²

で扱う計画だった。この形は悪い点が 3 種類(`z ∈ L`、`z ∈ L − w`、`z ≡ w`)あり、
分母が消える 2 種類は**それぞれ別の議論**を要する。

★★★本ブロックでは**分母のない形**——3 点 `P(z)`, `P(w)`, `−P(z+w)` の共線性の
**行列式**——に切り替える。

    collDet L w z
      = det [[℘(z), ℘'(z), 1], [℘(w), ℘'(w), 1], [℘(z+w), −℘'(z+w), 1]]
      = ℘'(w)(℘(z) − ℘(z+w)) − ℘(w)(℘'(z) + ℘'(z+w))
        + (℘(z)℘'(z+w) + ℘'(z)℘(z+w))

★これは `℘`・`℘'` の**積と和だけ**でできているので、悪い点は極の 2 種類
(`z ∈ L` と `z ∈ L − w`)しかない。さらに `z ↦ −z − w` は 3 点のうち 2 点を
入れ替えるので

    collDet L w (−z − w) = −collDet L w z            (`collDet_neg`)

であり、**第 2 種は第 1 種から反対称性だけで出る**。しかも Liouville で定数と分かった
あと、この反対称性から**その定数が 0 であること**まで出る(値を一点で計算しなくてよい)。

## ★★★★★★極の相殺は Taylor 3 項で `ring` に落ちる

`l₀ ∈ L` のまわりで `t := z − l₀`、`f(s) := ℘(s + w)` とおくと、周期性から
`℘(z + w) = f(t)`、`℘'(z + w) = f'(t)` であり、第 244 の `exists_pole_form` から

    ℘(z) = 1/t² + E(z),   ℘'(z) = −2/t³ + D(z)

なので、`collDet` の特異部は

    (2f(0) − 2f(t))/t³ + (f'(0) + f'(t))/t²
      = [ −2(f(t) − f(0)) + t(f'(t) + f'(0)) ] / t³

となる。★★★分子が `t³` で割り切れることが在庫の
`AnalyticAt.exists_eq_sum_add_pow_mul`(解析的な剰余つき Taylor 展開)で
**`ring` に落ちる**:`f` を 3 次まで、`f'` を 2 次まで展開して代入すると
1 次と 2 次の項がちょうど消える(`exists_taylor_pole_cancel`)。

★`AnalyticAt.order` を持ち出して「3 位で消える」と言う必要はない。
第 244 で採った「位数を数えず因数分解で」の方針がそのまま通る。

## ★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_taylor_pole_cancel` | ★★★★★★**Taylor 3 項で極が消える** |
| `shiftP`・`deriv_shiftP` 他 | ★平行移動した `℘` とその微分 |
| `collDet` | ★★★★★★**共線性の行列式** |
| `collDet_neg` | ★★★★反対称性(`z ↦ −z − w`) |
| `collDet_add_lattice` | ★★周期性 |
| `exists_analytic_collDet` | ★★★★★★★**格子点で極が相殺する** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real PeriodPair

/-! ## ★★★★★★Taylor 3 項で極が消える -/

/-- ★★★★★★**Taylor 3 項で極が消える**——`−2(f z − f 0) + z(f' z + f' 0)` は `z³` で
割り切れる。

`f` を 3 次まで、`f'` を 2 次まで解析的な剰余つきで展開すると、0 次・1 次・2 次の項が
ちょうど打ち消し合う。★**位数を数えずに `ring` で出る**のがこの形の利点である。 -/
theorem exists_taylor_pole_cancel (f : ℂ → ℂ) (hf : AnalyticAt ℂ f 0)
    (hf' : AnalyticAt ℂ (deriv f) 0) :
    ∃ H : ℂ → ℂ, AnalyticAt ℂ H 0 ∧
      ∀ z : ℂ, -2 * (f z - f 0) + z * (deriv f z + deriv f 0) = z ^ 3 * H z := by
  obtain ⟨F, hF, hFeq⟩ := hf.exists_eq_sum_add_pow_mul 3
  obtain ⟨G, hG, hGeq⟩ := hf'.exists_eq_sum_add_pow_mul 2
  refine ⟨fun z => G z - 2 * F z, hG.sub (analyticAt_const.mul hF), ?_⟩
  intro z
  have h1 := hFeq z
  have h2 := hGeq z
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, iteratedDeriv_zero,
    iteratedDeriv_one, smul_eq_mul, zero_add, pow_zero, pow_one, Nat.factorial] at h1 h2
  rw [show (2 : ℕ) = 1 + 1 from rfl, iteratedDeriv_succ'] at h1
  simp only [iteratedDeriv_one] at h1
  rw [h1, h2]
  push_cast
  ring

/-! ## ★平行移動した `℘` -/

/-- 平行移動した `℘`——`s ↦ ℘(s + w)`。 -/
noncomputable def shiftP (L : PeriodPair) (w : ℂ) : ℂ → ℂ := fun s => L.weierstrassP (s + w)

theorem shiftP_apply (L : PeriodPair) (w s : ℂ) : shiftP L w s = L.weierstrassP (s + w) := rfl

theorem deriv_shiftP (L : PeriodPair) (w : ℂ) :
    deriv (shiftP L w) = fun s => L.derivWeierstrassP (s + w) := by
  funext s
  show deriv (fun u : ℂ => L.weierstrassP (u + w)) s = _
  rw [deriv_comp_add_const, L.deriv_weierstrassP]

theorem analyticAt_shiftP (L : PeriodPair) (w s : ℂ) (h : s + w ∉ L.lattice) :
    AnalyticAt ℂ (shiftP L w) s :=
  AnalyticAt.comp (f := fun u : ℂ => u + w) (x := s)
    (L.analyticOnNhd_weierstrassP (s + w) h) (analyticAt_id.add analyticAt_const)

theorem analyticAt_deriv_shiftP (L : PeriodPair) (w s : ℂ) (h : s + w ∉ L.lattice) :
    AnalyticAt ℂ (deriv (shiftP L w)) s := by
  rw [deriv_shiftP]
  exact AnalyticAt.comp (f := fun u : ℂ => u + w) (x := s)
    (L.analyticOnNhd_derivWeierstrassP (s + w) h) (analyticAt_id.add analyticAt_const)

/-! ## ★★★★★★共線性の行列式 -/

/-- ★★★★★★**共線性の行列式**——3 点 `P(z)`, `P(w)`, `−P(z+w)` が一直線上にある条件。

    det [[℘(z), ℘'(z), 1], [℘(w), ℘'(w), 1], [℘(z+w), −℘'(z+w), 1]]

★分母がないので、悪い点は極の 2 種類(`z ∈ L` と `z ∈ L − w`)だけである。 -/
noncomputable def collDet (L : PeriodPair) (w z : ℂ) : ℂ :=
  L.derivWeierstrassP w * (L.weierstrassP z - L.weierstrassP (z + w))
    - L.weierstrassP w * (L.derivWeierstrassP z + L.derivWeierstrassP (z + w))
    + (L.weierstrassP z * L.derivWeierstrassP (z + w)
        + L.derivWeierstrassP z * L.weierstrassP (z + w))

/-- ★★★★**反対称性**——`z ↦ −z − w` は 3 点のうち 2 点を入れ替えるので符号が変わる。

★これが「第 2 種の悪い点」と「定数が 0 であること」を同時に片付ける。 -/
theorem collDet_neg (L : PeriodPair) (w z : ℂ) : collDet L w (-z - w) = -collDet L w z := by
  simp only [collDet]
  rw [show -z - w + w = -z by ring, show -z - w = -(z + w) by ring]
  simp only [L.weierstrassP_neg, L.derivWeierstrassP_neg]
  ring

/-- ★★**周期性**——`℘`・`℘'` が(ゴミ値も込みで)周期的なので全 `z` で成り立つ。 -/
theorem collDet_add_lattice (L : PeriodPair) (w : ℂ) (l : L.lattice) (z : ℂ) :
    collDet L w (z + (l : ℂ)) = collDet L w z := by
  simp only [collDet]
  rw [show z + (l : ℂ) + w = z + w + (l : ℂ) by ring]
  simp only [L.weierstrassP_add_coe, L.derivWeierstrassP_add_coe]

/-! ## ★★★★★★★格子点で極が相殺する -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★**格子点で極が相殺する**——`collDet L w` は `l₀ ∈ L` のまわりで
解析的な関数に一致する(`l₀` そのものを除く)。

`t := z − l₀`、`f := shiftP L w` として、特異部は

    (2f(0) − 2f(t))/t³ + (f'(0) + f'(t))/t²
      = [ −2(f(t) − f(0)) + t(f'(t) + f'(0)) ] / t³ = H(t)

となり、`exists_taylor_pole_cancel` で分子が `t³` で割り切れる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_analytic_collDet (L : PeriodPair) (w : ℂ) (hw : w ∉ L.lattice) (l₀ : L.lattice) :
    ∃ g : ℂ → ℂ, AnalyticAt ℂ g (l₀ : ℂ) ∧
      ∀ z : ℂ, z ≠ (l₀ : ℂ) → collDet L w z = g z := by
  obtain ⟨E, D, hE, hD, hP, hP'⟩ := exists_pole_form L l₀
  have hw0 : (0 : ℂ) + w ∉ L.lattice := by simpa using hw
  have hf0 : AnalyticAt ℂ (shiftP L w) 0 := analyticAt_shiftP L w 0 hw0
  have hf0' : AnalyticAt ℂ (deriv (shiftP L w)) 0 := analyticAt_deriv_shiftP L w 0 hw0
  obtain ⟨H, hH, hHeq⟩ := exists_taylor_pole_cancel (shiftP L w) hf0 hf0'
  have hQ : ∀ z : ℂ, L.weierstrassP (z + w) = shiftP L w (z - (l₀ : ℂ)) := by
    intro z
    rw [shiftP_apply, show z - (l₀ : ℂ) + w = z + w + ((-l₀ : L.lattice) : ℂ) by push_cast; ring,
      L.weierstrassP_add_coe]
  have hQ' : ∀ z : ℂ, L.derivWeierstrassP (z + w) = deriv (shiftP L w) (z - (l₀ : ℂ)) := by
    intro z
    simp only [deriv_shiftP]
    rw [show z - (l₀ : ℂ) + w = z + w + ((-l₀ : L.lattice) : ℂ) by push_cast; ring,
      L.derivWeierstrassP_add_coe]
  have hc : shiftP L w 0 = L.weierstrassP w := by rw [shiftP_apply, zero_add]
  have hc' : deriv (shiftP L w) 0 = L.derivWeierstrassP w := by
    simp only [deriv_shiftP, zero_add]
  have hshift : AnalyticAt ℂ (fun z : ℂ => z - (l₀ : ℂ)) (l₀ : ℂ) :=
    analyticAt_id.sub analyticAt_const
  have hH0 : AnalyticAt ℂ H ((l₀ : ℂ) - (l₀ : ℂ)) := by rw [sub_self]; exact hH
  have hf00 : AnalyticAt ℂ (shiftP L w) ((l₀ : ℂ) - (l₀ : ℂ)) := by rw [sub_self]; exact hf0
  have hf00' : AnalyticAt ℂ (deriv (shiftP L w)) ((l₀ : ℂ) - (l₀ : ℂ)) := by
    rw [sub_self]; exact hf0'
  have hHl : AnalyticAt ℂ (fun z : ℂ => H (z - (l₀ : ℂ))) (l₀ : ℂ) :=
    AnalyticAt.comp (f := fun z : ℂ => z - (l₀ : ℂ)) (x := (l₀ : ℂ)) hH0 hshift
  have hfl : AnalyticAt ℂ (fun z : ℂ => shiftP L w (z - (l₀ : ℂ))) (l₀ : ℂ) :=
    AnalyticAt.comp (f := fun z : ℂ => z - (l₀ : ℂ)) (x := (l₀ : ℂ)) hf00 hshift
  have hfl' : AnalyticAt ℂ (fun z : ℂ => deriv (shiftP L w) (z - (l₀ : ℂ))) (l₀ : ℂ) :=
    AnalyticAt.comp (f := fun z : ℂ => z - (l₀ : ℂ)) (x := (l₀ : ℂ)) hf00' hshift
  refine ⟨fun z => H (z - (l₀ : ℂ))
      + L.derivWeierstrassP w * (E z - shiftP L w (z - (l₀ : ℂ)))
      - L.weierstrassP w * (D z + deriv (shiftP L w) (z - (l₀ : ℂ)))
      + (E z * deriv (shiftP L w) (z - (l₀ : ℂ)) + D z * shiftP L w (z - (l₀ : ℂ))), ?_, ?_⟩
  · exact ((hHl.add (analyticAt_const.mul (hE.sub hfl))).sub
      (analyticAt_const.mul (hD.add hfl'))).add ((hE.mul hfl').add (hD.mul hfl))
  · intro z hz
    have ht : z - (l₀ : ℂ) ≠ 0 := sub_ne_zero.2 hz
    have hkey := hHeq (z - (l₀ : ℂ))
    rw [hc, hc'] at hkey
    simp only [collDet]
    rw [hP z, hP' z, hQ z, hQ' z]
    field_simp
    linear_combination hkey

/-! ## ★出典の紐付け(`.src`) -/

def exists_taylor_pole_cancel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——Taylor 3 項で極が消える)",
    sectionId := "genell-def-3-3" }

def collDet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——共線性の行列式)",
    sectionId := "genell-def-3-3" }

def exists_analytic_collDet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——格子点で極が相殺する)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
