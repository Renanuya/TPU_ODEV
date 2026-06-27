# TPU_ODEV — Ortak Emetörlü Yükselteç Simülatörü

Proje 4: Ortak emetörlü (OE) yükselteç devresinin **DC analizi** ve **AC küçük işaret (hibrit-π parametre)** analizini gerçekleştiren web tabanlı simülatör.

**Canlı demo:** https://site-sage-pi-25.vercel.app

## Özellikler

Arayüz dikey üç bloktan oluşur:

1. **Ortak Emetörlü Yükselteç Devresi** — şema + tüm eleman değerleri kaydırıcı (araç çubuğu): `Vcc, R1, R2, RC, RE, RL, Rs, β, vs, frekans`.
2. **DC Eş Değer Devre** — şema + hesaplanan `IB, IC (ICQ), IE, VCE, VCB, VBE` ve çalışma bölgesi (**AKTİF / DOYMA / KESİM**).
3. **AC Eş Değer Devre (Hibrit-π)** — şema + `Kv, Kv0, Ki, Kvg, Rg (giriş direnci), Ro, gm, rπ` ve giriş–çıkış gerilim grafiği (`vi(t)` & `vo(t)` üst üste, 180° ters faz, kırpılma uyarısı).

## Teori

- DC: gerilim böleni Thevenin → `IB = (VTh − VBE) / [RTh + (β+1)·RE]`, `IC = β·IB`, `VCE = VCC − IC·RC − IE·RE`.
- Bölge kararı: `VTh ≤ 0.5 V → KESİM`; `β·IB ≥ IC(sat) → DOYMA`; aksi halde AKTİF.
- AC (orta bant, ro→∞): `gm = ICQ/VT`, `rπ = β/gm`, `Rg = R1∥R2∥rπ`, `Kv = −gm·(RC∥RL)`, `Kvs = Kv·Rg/(Rg+Rs)`.

Formüller raporun sayısal doğrulama örneğiyle birebir tutar (Kv=−91.8, Kvs=−71.9 vb.) ve 4 farklı devre (AKTİF/DOYMA/KESİM) ile çapraz kontrol edilmiştir.

## Çalıştırma

Saf statik site, build yok:

```bash
cd site
python -m http.server 8077
# http://localhost:8077
```

## Yapı

```
site/        statik simülatör (index.html, styles.css, app.js)
docs/        proje destek raporu (PDF)
```

Kaynak: Sedra/Smith *Microelectronic Circuits*, Boylestad *Electronic Devices*, IIT Bombay OE amplifier ders notu.
